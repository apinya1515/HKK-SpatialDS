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

library(nimble) ## this on code can call nimble package
library(dplyr)
library(coda) # for post-training mcmc check 
library(MCMCvis)
library(sf) # for import spatial data
library(spdep) # for defining spatial neighbors (used for CAR)

# Transect data (table 1)
#data_tr <- read.table('D:\\UserData\\Dropbox\\@KU_pop_course\\day1_exercise\\line_data.txt', 
#                       sep='\t', header=T)
data_tr <- read.csv('D:\\2024SpatialDistance\\analysis\\hkkdata.csv')#Apinya
# landscape data
#data_land <- read.csv('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\HKK_Cov1sqkm_.csv')
data_land <- read.csv('D:\\2024SpatialDistance\\analysis\\HKK_Cov1sqkm_15042025.csv')#Apinya
# transect X landscape data
#data_prop <- read.csv('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\TRidentity.csv')
data_prop <- read.csv('D:\\2024SpatialDistance\\analysis\\TRidentity.csv')#Apinya
##############################################################################
########################## INITIALIZATION ####################################
### Input parameters
species <- 'SBR'
dist_limit <- 120 # max observed distance
dist_class_n <- 3 # number of dist classes
gsBreaks <- c(1, 2, 3) # upper breaks group size of each gs classes - c(1, 3) = 2 classes 1st = 1, 2nd >= 2
##############################################################################

### Get data from files
tran_id <- sort(unique(data_tr$Tr.no)) # name of transect
tran_n <- length(tran_id) # number of transect
grid_id <- data_land$grid_id # grid id
grid_n <- nrow(data_land) # number of grid
covar <- as.data.frame(data_land[,-1] %>% scale()) # select only covariate columns and rescale

### Create observation matrix (matrix=transect, row=cluster size class ,column=distance class)
dist_width <- dist_limit/dist_class_n # interval between distance class
distBreaks = seq(dist_width, dist_class_n*dist_width, by=dist_width) # upper breaks of distance classes

# subset data by species
data_sub <- data_tr %>% filter(Species==species)

gs_max <- max(data_sub$Gz.sz)   # max group size
gs_class_n <- length(gsBreaks) # number og gs classes

#### Create zeros-array for the data
y_matrix <- array(0, dim=c(gs_class_n, dist_class_n, tran_n), 
                  dimnames=list(1:gs_class_n, distBreaks, tran_id)) # 3 dimensions = group size, distance, transect

### Populate y data for target species ... by +1 to y_matrix for each data point that fall into the cell
# y_matrix[group_size, distance, transect]
for(i in 1:nrow(data_sub)){
  samp <- data_sub[i,] # get data
  if(samp$P.dist > dist_limit){print(paste0('data ',i, ' distance:',  round(samp$P.dist,2), ' OUT OF RANGE !!!'));next;}
  tr <- which(tran_id ==  samp$Tr.no) # transect index
  j <- sum(distBreaks <= samp$P.dist) + 1 # distance class index
  k <- min(which(gsBreaks >= samp$Gz.sz)) # group size class index
  print(paste0('data ',i, ' gs:', samp$Gz.sz, ' distance:',  round(samp$P.dist,2),' in grsz:', k, ' dist:', j, ' transect:', tr))
  y_matrix[k, j, tr] = y_matrix[k, j, tr] + 1 # add observation into group class - dist class - transect
}
print(y_matrix)

sum(y_matrix) # check the number of observations
nrow(data_sub) # check the number of observations

### Grid - Transect Mapping ###
# Create a matrix row=transect col=grid to map proportion from each grid to transect
head(data_prop)
lenM <- matrix(0, nrow=tran_n, ncol=grid_n) # create zeo matrix row=#of trasect col=#of grid
for(i in 1:nrow(data_prop)){
  dat <- data_prop[i,]
  lenM[dat$tr_seq, dat$grid_id] <- dat$Length
}
# total length for each transect
lenSum <- rowSums(lenM)
# create proportion matrix
propM <- sweep(lenM, 1, lenSum, FUN="/")
propM[1:30, 1:30]
rowSums(propM) # check that each row sum to 1

### Spatial dependency ######################
### CAR elements #############################
#poly <- st_read('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\shp\\HKK1sqkmGrid.shp')
poly <- st_read('D:\\2024SpatialDistance\\analysis\\HKK1sqkmGrid.shp')#Apinya
poly$ID <- seq(1:grid_n)
# Create a neighbors list (Queen's contiguity)
nb <- poly2nb(poly, queen = TRUE, row.names = poly$ID)
adj <- unlist(nb) # adjacency index (grid id of each neighbor)
njoin <- length(adj) # number of spatial joins
num <- sapply(nb, length) # number of neighbors of each grid 
### (sum of num must match with the length of adj)

####################################################################################
####################################################################################
####################################################################################
# Model code in NIMBLE
distanceModelCode <- nimbleCode({
  ##### Priors #####
  for(i in 1:n_covar){
    beta[i] ~ dnorm(0, 0.001) # priors of landscape var regression coefficient
    w[i] ~ dbern(0.5) # priors of binary indicator for each coefficient
  }
  muc ~ dunif(1, 3) # average cluster size
  #sigma0 ~ dunif(1, 10) # sigma0 # try larger prior to avoid -Inf logProb error
  sigma0 ~ dlnorm(2, 0.5)  # Mean ~7.4, allows larger values
  p ~ dunif(0.1, 10) # model parameter for sigma ~ group_size
  
  # CAR priors
  tau ~ dgamma(0.001, 0.001) # precision param for spatial CAR
  weights[1:njoin] <- 1 # weight of all spatial joins = 1
  b_spatial[1:L] ~ dcar_normal(adj = adj[1:njoin], weights = weights[1:njoin], 
                           num = num[1:L], tau = tau, zero_mean = 0) # random spatial effect (intercept) intercept for each grid
  
  ##### Model #####
  
  ### Landscape abundance for each grid
  w_beta[1:n_covar] <- w[1:n_covar] * beta[1:n_covar]
  for(l in 1:L) {
    # Regression abundance + b_spatial(CAR elements)
    log_lam[l] <- inprod(covar[l, 1:n_covar], w_beta[1:n_covar]) + b_spatial[l] # regression for log abundance (without intercept)
    lam[l] <- exp(log_lam[l]) # lambda for each grid
  }
  
  ### Abundance at each transect
  for (i in 1:I) {
    tr_lam[i] <- inprod(propM[i,1:L], lam[1:L])
  }
  
  ### Model prob group size class
  for(m in 1:gs_max) { # loop through every group size
    ### calculate truncated Poisson (at m=1) prob for each size class from 1 to gs_max
    gs_m[m] <- (exp(-muc) * pow(muc, m)) / (factorial(m) * (1 - exp(-muc)))
  }
  ### calculate prob of each group size class & mean group size for each class
  if(gsBreaks[1] == 1) { # if the first class contain only group size (m) = 1
    gs_k[1] <- gs_m[1]  # Direct assignment for single element
    gs_k_mean[1] <- gs_m[1] * 1 / gs_k[1]  # Mean for class 1
  } else {
    gs_k[1] <- sum(gs_m[1:gsBreaks[1]])  # Sum for multiple elements
    gs_k_mean[1] <- sum(gs_m[1:gsBreaks[1]] * 1:gsBreaks[1]) / gs_k[1]  # Mean for class 1 (see page 45)
  }
  for(k in 2:K) { # loop through group size class 2 and on
    gs_k[k] <- sum(gs_m[(gsBreaks[k-1]+1):gsBreaks[k]]) # sum of gs_m for class 2 and on
    gs_k_mean[k] <- sum(gs_m[(gsBreaks[k-1]+1):gsBreaks[k]] * (gsBreaks[k-1]+1):gsBreaks[k]) /  gs_k[k] # gs mean for class 2 and on (see page 45)
  }
  
  ### Calculate sigma for each size class
  # sigma is function of group size class (k) bigger group size -> larger sigma[k] -> slower decay half normal detection function
  for(k in 1:K){
    log(sigma[k]) <- sigma0 + p * (gs_k_mean[k] - 1)
  }
  
  ### Distance sampling model
  # Compute pi for each distance class [j] and group size class [k]
  # pi[k.j] = probability of detection is depend on group size(k) and distance(j) but is the same across transect
  for(k in 1:K){
    # Compute the CDF at dist_limit for normalization
    # Normalization: Added F_dist_limit[k] as the CDF at dist_limit (i.e., distBreaks[J]), and normalized mn_cell[k, j] by dividing by this value.
    # Used pnorm (NIMBLE’s standard normal CDF) instead of phi for clarity and compatibility.
    # 2 * pnorm to uses a half-normal detection function for distance sampling
    F_dist_limit[k] <- 2 * pnorm(distBreaks[J] / sigma[k]) - 1
    # Distance class 1
    mn_cell[k, 1] <- (2 * pnorm(distBreaks[1] / sigma[k]) - 1) / F_dist_limit[k]
    pi[k, 1] <- mn_cell[k, 1] * gs_k[k]
    # Distance classes 2 to J
    for(j in 2:J){
      mn_cell[k, j] <- 2 * (pnorm(distBreaks[j] / sigma[k]) - pnorm(distBreaks[j - 1] / sigma[k])) / F_dist_limit[k]
      pi[k, j] <- mn_cell[k, j] * gs_k[k]
    }
  }

  #### Model data for each transect
  # Loop for each element in all data matrix
  for(i in 1:I){
    for(k in 1:K){
      for(j in 1:J){
        # expected observation in transect i, group size class k, and distance class j
        mu[k, j, i] <- tr_lam[i] * pi[k, j] # abundance of transect i * probability of class k and j
        ### Model transect data
        y[k, j, i] ~ dpois(mu[k, j, i])
      }
    }
  }
  
})

####################################################################################

constants <- list(I = tran_n, J = dist_class_n, 
                  K = gs_class_n, gs_max = gs_max,
                  L = grid_n, n_covar = ncol(covar),
                  distBreaks = distBreaks, gsBreaks = gsBreaks,
                  adj = adj, num = num, njoin = njoin)
data <- list(y = y_matrix, covar = covar, propM = propM)
inits <- list(muc = runif(1, 1, 3), sigma0 = runif(1, 2, 5), p = runif(1, 0.5, 1.5), tau = runif(1,0.5,1),
              beta = rnorm(ncol(covar), 0,1), b_spatial = runif(grid_n, 0, 0.1), 
              w = rep(1, ncol(covar)))

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, 
                             constants = constants, 
                             data = data, 
                             inits = inits)

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = c("sigma", "p", "gs_m","gs_k", "sigma0", "pi", 'lam', 'beta', 'tau', 'b_spatial', 'w'))
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

# Run the MCMC
t_start <- Sys.time() # start time
samples_chains <- runMCMC(Cmcmc, niter = 10000, nburnin = 5000, thin = 2, nchains = 2)
(t_elapse <- Sys.time() - t_start) # time elapse
###summary(samples_chains)

################################ Post training ############################
# Convergence check - use package coda
# Diagnostics
# Convert samples to coda format
mcmc.list <- as.mcmc.list(lapply(samples_chains, as.mcmc))
# Filter parameters with prefixes "beta", "sigma", "w", and "muc"
param_names <- colnames(mcmc.list[[1]])
selected_params <- param_names[grep("^(beta|sigma|w|muc)", param_names)]
colnames(covar) # see covariates names
# Subset mcmc.list to include only selected parameters
sub_mcmc.list <- mcmc.list[, selected_params, drop = FALSE]
MCMCtrace(sub_mcmc.list, pdf = FALSE) # Visualize MCMC trace
gelman.diag(sub_mcmc.list)  # R-hat (R-hat < 1.1 for convergence)
effectiveSize(sub_mcmc.list)  # ESS (ESS 100 - 200 for convergence)
autocorr.plot(sub_mcmc.list)  # Autocorrelation (Low Autocorrelation for convergence)
MCMCsummary(sub_mcmc.list)  # Summary with R-hat and ESS

# merge multiple chain into 1
samples <- do.call(rbind, samples_chains)

# fitted detection function
sigma_1 <- median(samples[,'sigma[1]'])
sigma_2 <- median(samples[,'sigma[2]'])
p <- median(samples[,'p'])
p

# multinomial pi
pi_val <- apply(samples[,grep("^pi", colnames(samples))], 2, median)
pi <- matrix(pi_val, nrow=gs_class_n, byrow = F)
pi
sum(pi)

### Visualize pi
plot(pi[1,], type='l', xlab='Distance Class', ylab='Group Size', main='Multinomial Prob')
lines(pi[2,], lty=2)
legend('topright', legend=c('1','2'), title = 'Group size', lty=c(1,2))

# calculate gx (only distance)
x_seq <- seq(0, dist_limit)
gx_1 <- exp((-1*x_seq^2)/(2*sigma_1^2))
gx_2 <- exp((-1*x_seq^2)/(2*sigma_2^2))
hp <- hist(data_sub$P.dist, main='Histogram & Detection function', xlab='Distance')
lines(x_seq, gx_1*max(hp$counts), col='red', lwd=2)
lines(x_seq, gx_2*max(hp$counts), col='green', lwd=2)
legend('topright', legend=c('1','2'), title = 'Group class', col=c('red', 'green'), lty=1)

# Site-specific Lambda (Abundance)
lam <- samples[,grep("^lam", colnames(samples))]
tr_names <- dimnames(y_matrix)[[3]]
colnames(lam) <- 1:grid_n
boxplot(lam[,1:50], ylim=c(0,100))
boxplot(lam[,101:150], ylim=c(0,100))
abun_med <- apply(lam, 2, median, na.rm = TRUE)
hist(abun_med)

# Site-specific random effect
b_spatial <- samples[,grep("^b_spatial", colnames(samples))]
bspat_med <- apply(b_spatial, 2, median, na.rm = TRUE)

# Beta
beta <- samples[,grep("^beta", colnames(samples))]
colnames(beta) <- colnames(covar)
boxplot(beta)
abline(h=0, lty=3)

# Indicator
w <- samples[,grep("^w", colnames(samples))]
colnames(w) <- colnames(covar)
w_contrib <- colSums(w)/nrow(w)
barplot(w_contrib)

# plot beta that w has value at least 50% of 1
beta_w <- beta[,w_contrib > 0.5]
boxplot(beta_w, ylim=c(min(0,min(beta_w)),max(0,max(beta_w))))
abline(h=0, lty=3)

# Inject back to shp file
library(terra)
abun_med[abun_med>20]
poly$abun_med <- abun_med
poly$bspat_med <- bspat_med
plot(poly["abun_med"], main='Median grid-level abundance (lambda)')
plot(poly["bspat_med"], main='CAR \n Median spatial random effects')

