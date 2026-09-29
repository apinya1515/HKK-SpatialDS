### Non-centred, sparse version of CAR-Kumar2021-IND.R ###

###### 29/09/2026
### Same statistical model as CAR-Kumar2021-IND.R (same priors, same likelihood); only the
### parameterisation and bookkeeping change, to fix slow mixing of sigma_spatial / b_spatial:
###   1. Non-centred CAR: u ~ dcar_normal(tau = 1, zero_mean = 1), b_spatial = sigma_spatial * u.
###      b_spatial ~ dcar_normal(tau = 1/sigma_spatial^2) exactly as before, but sigma_spatial is
###      no longer tied to 3013 mostly data-free random effects.
###   2. Sparse transect mapping: only the T grids crossed by a transect (propM > 0) enter the
###      likelihood, so z is computed for those grids only (zT). Previously every grid fed every
###      transect through inprod(propM[i, 1:L], z[1:L]), making each b_spatial update O(L * I).
###   3. w is a constant covariate mask (it was a fixed node with its sampler removed), and grid
###      abundance / TOTAL_ABUND are derived after sampling from beta0, beta, sigma_spatial, u, AGS.
### Requires (from @data_prepare_011025.R): tran_n, grid_n, covar, propM, water_mask, adj, num,
###   njoin, lenSum, logFactorial, dist/gs settings; plus w_mask (0/1 vector over covariates).

distanceModelCode <- nimbleCode({
  ##### Priors #####
  for (i in 1:n_covar) {
    beta[i] ~ dnorm(0, sd = 1.5) # weakly informative priors of landscape var regression coefficient
  }
  beta0 ~ dnorm(0, sd = 1.5) # weakly informative intercept prior

  muc ~ dunif(1, gs_max) # universal average cluster size (for truncated-Poisson)
  sigma0 ~ dnorm(3, sd = 2)
  p ~ dunif(0, 5) # sigma ~ group size (p >= 0: larger groups have larger detection range)

  # CAR priors (non-centred)
  sigma_spatial ~ dunif(0, 1.5) # SD of spatial CAR effect
  weights[1:njoin] <- 1
  u[1:L] ~ dcar_normal(
    adj = adj[1:njoin], weights = weights[1:njoin],
    num = num[1:L], tau = 1, zero_mean = 1
  ) # standardised spatial effect, sum-to-zero; b_spatial = sigma_spatial * u

  ##### Model #####

  ### Grid-level abundance, only for grids crossed by transects
  w_beta[1:n_covar] <- w[1:n_covar] * beta[1:n_covar]
  for (t in 1:T) {
    zT[t] <- exp(beta0 + inprod(covarT[t, 1:n_covar], w_beta[1:n_covar]) + sigma_spatial * u[gT[t]]) * maskT[t]
  }

  ### Abundance at each transect
  for (i in 1:I) {
    Z[i] <- inprod(propT[i, 1:T], zT[1:T]) # identical to inprod(propM[i, 1:L], z[1:L])
    area_tr[i] <- (2 * distBreaks[J] * lenSum[i]) / 1000000
    lam[i] <- nrep * Z[i] * area_tr[i]
  }

  ### Group size (zero-truncated Poisson)
  for (m in 1:gs_max) {
    summand[m] <- exp(m * log(muc) - logFactorial[m])
    summandForMean[m] <- m * summand[m]
  }
  C_ztp <- sum(summand[1:gs_max])

  if (gsBreaks[1] == 1) {
    gs_k[1] <- summand[1] / C_ztp
    gs_k_mean[1] <- summandForMean[1] / summand[1]
  } else {
    gs_k[1] <- sum(summand[1:gsBreaks[1]]) / C_ztp
    gs_k_mean[1] <- sum(summandForMean[1:gsBreaks[1]]) / sum(summand[1:gsBreaks[1]])
  }
  for (k in 2:K) {
    gs_k[k] <- sum(summand[(gsBreaks[k - 1] + 1):(gsBreaks[k])]) / C_ztp
    gs_k_mean[k] <- sum(summandForMean[(gsBreaks[k - 1] + 1):(gsBreaks[k])]) / sum(summand[(gsBreaks[k - 1] + 1):(gsBreaks[k])])
  }

  ### Detection scale for each group size class
  for (k in 1:K) {
    sigma[k] <- exp(sigma0 + p * (gs_k_mean[k] - 1))
  }

  ### Multinomial cell probabilities (group size class k x distance class j)
  for (k in 1:K) {
    mn_cell[k, 1] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * (pnorm(distBreaks[1], mean = 0, sd = sigma[k]) - 0.5)
    pi[k, 1] <- gs_k[k] * mn_cell[k, 1]
    for (j in 2:J) {
      mn_cell[k, j] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * (pnorm(distBreaks[j], mean = 0, sd = sigma[k]) - pnorm(distBreaks[j - 1], mean = 0, sd = sigma[k]))
      pi[k, j] <- gs_k[k] * mn_cell[k, j]
    }
  }

  #### Observations
  for (i in 1:I) {
    for (k in 1:K) {
      for (j in 1:J) {
        mu[k, j, i] <- lam[i] * pi[k, j]
        y[k, j, i] ~ dpois(mu[k, j, i])
      }
    }
  }

  ### Derived
  AGS <- muc / (1 - exp(-muc))
})

# Grids crossed by at least one transect
gT <- which(colSums(propM) > 0)

tracked_var <- c("beta0", "beta", "muc", "sigma0", "sigma", "p", "pi", "gs_k", "AGS", "sigma_spatial", "u")

constants <- list(
  nrep = nrep, I = tran_n, J = dist_class_n,
  K = gs_class_n, gs_max = gs_max,
  L = grid_n, T = length(gT), gT = gT, n_covar = ncol(covar),
  distBreaks = distBreaks, gsBreaks = gsBreaks,
  adj = adj, num = num, njoin = njoin,
  lenSum = lenSum, logFactorial = logFactorial
)
data <- list(
  y = y_matrix, w = w_mask,
  covarT = as.matrix(covar)[gT, , drop = FALSE],
  propT = propM[, gT, drop = FALSE],
  maskT = water_mask[gT]
)
# sigma0 starts >= 3 (detection scale >= 20 m): with sigma0 near 2 the far distance-class
# probabilities underflow to 0, the initial logProb is -Inf, and the slice sampler on
# sigma_spatial then accepts any value (it left its prior support in BTG Rank 2, chain 1).
inits <- list(
  muc = runif(1, 1, gs_max), sigma0 = runif(1, 3, 5), p = runif(1, 0, 0.1),
  beta0 = 0, beta = rnorm(ncol(covar), 0, 0.5),
  u = { u0 <- rnorm(grid_n, 0, 0.1); u0 - mean(u0) },
  sigma_spatial = runif(1, 0.3, 1.2)
)

# Grid-level individual abundance for all L grids from a matrix of samples (post hoc).
# Returns list(b_spatial = draws x L, TOTAL_ABUND = draws).
derive_abundance <- function(smp, covar_all, w_mask, water_mask) {
  beta_cols <- sprintf("beta[%d]", seq_len(ncol(covar_all)))
  u_cols <- sprintf("u[%d]", seq_len(nrow(covar_all)))
  w_beta <- sweep(smp[, beta_cols, drop = FALSE], 2, w_mask, "*")        # draws x n_covar
  b_sp <- smp[, "sigma_spatial"] * smp[, u_cols, drop = FALSE]           # draws x L
  lin <- tcrossprod(w_beta, as.matrix(covar_all)) + smp[, "beta0"] + b_sp  # draws x L
  z <- sweep(exp(lin), 2, water_mask, "*")
  list(b_spatial = b_sp, TOTAL_ABUND = rowSums(z) * smp[, "AGS"])
}
