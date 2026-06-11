### Non-indicator version ###

###### Based on Kumar, 2021
###### 12/04/2025 ######
### Add binary indicator ("w") w/ Bernoulli for each regression coefficient
### Add group size classes (defined by "gsBreaks" in INITIALIZATION)
### ## calculate prob of each groups size ("gs_m") and prob of each group size class ("gs_k")
### ## "gs_k" is used to calculate "pi" (product between "gs_k" and cond.prob of each distance class ("mn_cell[k, j]"))
### ## expected size of each group size class ("gs_k_mean") is used to calculate "sigma"
### Add post-training MCMC diagnostics; R-hat, ESS, and autocorrelation (for nchains >= 2 only)
###### 22/04/2025 ######
### Change custom half-normal detection prob to Nimble-compatibility function pnorm (with 2 multiplication)
### Calculate cumulative prob at largest ditance (F_dist_limit[k]) and normalized other distance class with this
### Change prior of sigma0 to lognormal to allow larger values

###### 22/07/2025
### Split code into 2 files:
###   main_run_model = main file for pre-model and post-model
###   model*** = nimble model code to be called by source() function
### Changes:
###   Larger sigma0 prior to avoid warning of too small pi
###   Changed function of zero-truncated Poisson to be the same as in the book
###   Changed half-normal function to be the same as in the book
###   Add grid_int = grid-level non-spatial intercept for cluster abundance

###### 30/09/2025
### Separate fixed-effect (fix_lam) and spatial random effect for grid-level abundance to made it easier to be tracked
### models will have non-indicator version

###### 11/11/2025
### Add site(transect-level) specific abundance
### remove non-spatial intercept from grid-level abundance model
### Change grid-level abundance from 'lam' to 'z'
### Change transect-level abundance from 'tr_lam' to 'lam'
### Add transect-level intercept 'alpha' to function of log(lam)
### Change ABUND formula to: "ABUND[l] <- z[l] * AGS"

####################################################################################
####################################################################################
####################################################################################
# Model code in NIMBLE
distanceModelCode <- nimbleCode({
  ##### Priors #####
  for (i in 1:n_covar) { # Loop through covariates to create priors for coefficients (beta) and indicator (w)
    beta[i] ~ dnorm(0, sd = 1.5) # weakly informative priors of landscape var regression coefficient
  }
  # Site(transect)-specific effect on abundance
  alpha ~ dnorm(0, sd = 0.5) # weakly informative prior to avoid confounding with intercept

  muc ~ dunif(1, gs_max) # universal average cluster size (for truncated-Poisson)
  # sigma0 ~ dunif(1, 10) # sigma0 # try larger prior to avoid -Inf logProb error
  sigma0 ~ dnorm(3, sd = 2) # ***larger values - to prevent too low logProb (less than -1e12)
  p ~ dnorm(0, sd = 1) # model parameter for sigma ~ group_size

  # CAR priors
  sigma_spatial ~ dunif(0, 1.5) # weakly informative standard deviation prior for spatial CAR
  tau <- 1 / (sigma_spatial^2) # precision param for spatial CAR
  weights[1:njoin] <- 1 # weight of all spatial joins = 1
  b_spatial[1:L] ~ dcar_normal(
    adj = adj[1:njoin], weights = weights[1:njoin],
    num = num[1:L], tau = tau, zero_mean = 0
  ) # random spatial effect (intercept) intercept for each grid

  ##### Model #####

  ### Landscape abundance for each grid
  for (l in 1:L) { # Loop trough each (of all) grid
    # Regression for mean abundance with b_spatial(CAR elements)
    # z = grid-level abundance
    fix_z[l] <- inprod(covar[l, 1:n_covar], beta[1:n_covar]) # regression for grid-level mean abundance
    z[l] <- exp(fix_z[l] + b_spatial[l]) * water_mask[l] # mean abundance as sum of fixed effect and spatial random effect, masked for water
  }

  ### Abundance at each transect
  # lam = transect-level abundance
  for (i in 1:I) { # Loop through each transect
    # Z[i] is proportionated density along transect (groups / sq km)
    Z[i] <- inprod(propM[i, 1:L], z[1:L])
    # Search area of the transect strip in square km (since grid cells are 1 sq km)
    # total width is 2 * distBreaks[J] meters, length is lenSum[i] meters
    area_tr[i] <- (2 * distBreaks[J] * lenSum[i]) / 1000000
    # Expected number of groups in strip: repeat walks * density * search area * site alpha
    lam[i] <- nrep * Z[i] * area_tr[i] * exp(alpha)
  }


  ### Modeling prob for each group size cateogory
  for (m in 1:gs_max) { # loop through every group size
    # !!!## Array storing probability mass function (PMF) values for the Poisson distribution.
    # !!!## summand[m] = exp(m * log(muc) - logFactorial[m]) = exp(log(muc^m) - log(m!)) = (muc^m)/m!
    summand[m] <- exp(m * log(muc) - logFactorial[m])
    # !!!## Array storing terms used to compute the expected value of the zero-truncated Poisson.
    # !!!## summandForMean[m] = m*summand[m] = m*(muc^m)/m! = (muc^m)/(m-1)!
    summandForMean[m] <- m * summand[m]
  }
  C_ztp <- sum(summand[1:gs_max]) # The sum of unnormalized probabilities

  ### calculate prob of each group size class & mean group size for each class
  # gs_k[k] = probability that animal cluster is size category k
  # gs_k_mean[k] = mean cluster size of size category k
  if (gsBreaks[1] == 1) { # if first distance category breaks is 1
    gs_k[1] <- summand[1] / C_ztp # prob for category 1
    gs_k_mean[1] <- summandForMean[1] / summand[1]
  } else { # else if first distance category breaks is >= 2
    gs_k[1] <- sum(summand[1:gsBreaks[1]]) / C_ztp # prob for category 1
    gs_k_mean[1] <- sum(summandForMean[1:gsBreaks[1]]) / sum(summand[1:gsBreaks[1]])
  }
  for (k in 2:K) { # prob for group category >= 2
    gs_k[k] <- sum(summand[(gsBreaks[k - 1] + 1):(gsBreaks[k])]) / C_ztp
    gs_k_mean[k] <- sum(summandForMean[(gsBreaks[k - 1] + 1):(gsBreaks[k])]) / sum(summand[(gsBreaks[k - 1] + 1):(gsBreaks[k])])
  }

  ### *****
  ### Calculate sigma for each size class
  # sigma is function of group size class (k) bigger group size -> larger sigma[k] -> slower decay half normal detection function
  for (k in 1:K) { # loop through group size categories
    sigma[k] <- exp(sigma0 + p * (gs_k_mean[k] - 1))
  }

  ### Distance sampling model
  # Compute pi for each distance class [j] and group size class [k]
  # pi[k.j] = probability of detection is depend on group size(k) and distance(j) but is the same across transect
  for (k in 1:K) {
    # Compute half-normal detection function based on sigma[k] for each cluster size category k
    # The calculation is based on CDF of normal distribution using function "phi()" /// phi(y) if y=x/sd then it is standard normal
    # multiplying with the first term "sqrt(2*pi)*sigma[k]" to undoes the scaling of standard normal distribution, and adjust to match with sigma[k]
    # The dividing by "distBreaks[K]" is to normalize the prob of based on maximum distance
    # Distance class 1
    # Distance class 1: normalize by distBreaks[J] (maximum distance class J)
    mn_cell[k, 1] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * (pnorm(distBreaks[1], mean = 0, sd = sigma[k]) - 0.5)
    # calculate multinomial detection prob (gs[k] * gs[k,j])
    pi[k, 1] <- gs_k[k] * mn_cell[k, 1]
    # Distance classes 2 to J (through loop)
    for (j in 2:J) {
      # calculate difference between CDF at break of j class and j-1 class, normalize by distBreaks[J]
      mn_cell[k, j] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * (pnorm(distBreaks[j], mean = 0, sd = sigma[k]) - pnorm(distBreaks[j - 1], mean = 0, sd = sigma[k]))
      # calculate multinomial detection prob (gs[k] * gs[k,j])
      pi[k, j] <- gs_k[k] * mn_cell[k, j]
    }
  }

  #### Model data for each transect
  # Loop for each element in all data matrix
  for (i in 1:I) {
    for (k in 1:K) {
      for (j in 1:J) {
        # expected observation in transect i, group size class k, and distance class j
        mu[k, j, i] <- lam[i] * pi[k, j] # abundance of transect i * probability of class k and j
        ### Model transect data
        y[k, j, i] ~ dpois(mu[k, j, i])
      }
    }
  }

  ##########################
  ### Derived parameters ###
  # average group size
  AGS <- muc / (1 - exp(-muc))
  # abundance for each grid
  ### ABUND[1:L] <- z[1:L] * AGS
  for (l in 1:L) {
    ABUND[l] <- z[l] * AGS # generate abundance for each grid from Poisson distn with z
  }
  # Total abundance
  TOTAL_ABUND <- sum(ABUND[1:L])
})

# Tracked variables' names
tracked_var <- c(
  "sigma", "p", "muc", "gs_k", "sigma0", "pi", "fix_z", "z", "Z", "lam",
  "beta", "alpha", "tau", "fix_z", "b_spatial", "AGS", "ABUND",
  "TOTAL_ABUND"
)

constants <- list(
  nrep = nrep, I = tran_n, J = dist_class_n,
  K = gs_class_n, gs_max = gs_max,
  L = grid_n, n_covar = ncol(covar),
  distBreaks = distBreaks, gsBreaks = gsBreaks,
  adj = adj, num = num, njoin = njoin,
  lenSum = lenSum, logFactorial = logFactorial
)
data <- list(y = y_matrix, covar = covar, propM = propM, water_mask = water_mask)
inits <- list(
  muc = runif(1, 1, gs_max), sigma0 = runif(1, 2, 5), p = 0,
  beta = rnorm(ncol(covar), 1, 2), b_spatial = runif(grid_n, 0, 0.1),
  alpha = runif(1, -0.5, 0.5), sigma_spatial = runif(1, 0.5, 1.2)
)
# saveRDS(constants, 'constants.RDS')
# saveRDS(data, 'data.RDS')
# saveRDS(inits, 'inits.RDS')
