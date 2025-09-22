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

library(nimble)
library(dplyr)
library(coda) # for post-training mcmc check 
library(MCMCvis)
library(sf) # for import spatial data
library(spdep) # for defining spatial neighbors (used for CAR)

# Transect data (table 1)
data_tr <- read.table('D:\\UserData\\Dropbox\\@KU_pop_course\\day1_exercise\\line_data.txt', 
                       sep='\t', header=T)
# landscape data
data_land <- read.csv('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\HKK_Cov1sqkm_.csv')
# transect X landscape data
data_prop <- read.csv('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\TRidentity.csv')

##############################################################################
########################## INITIALIZATION ####################################
### Input parameters
species <- 'SBR'
dist_limit <- 120 # max observed distance
data_tr$P.dist[data_tr$P.dist > dist_limit ] <- dist_limit  # set maximum distance to the limit
gsBreaks <- c(1, 2, 3) # upper breaks group size of each gs classes - c(1, 3) -> 2 classes 1st class = 1, 2nd class >= 2, the maximum group size is 3
dist_class_n <- 5 # number of dist classes
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
  tr <- which(tran_id ==  samp$Tr.no) # transect index
  j <- sum(distBreaks < samp$P.dist) + 1 # distance class index
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
poly <- st_read('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\shp\\HKK1sqkmGrid.shp')
poly$ID <- seq(1:grid_n)
# Create a neighbors list (Queen's contiguity = 8 neighbors)
nb <- poly2nb(poly, queen = TRUE, row.names = poly$ID)
adj <- unlist(nb) # adjacency index (grid id of each neighbor)
njoin <- length(adj) # number of spatial joins
num <- sapply(nb, length) # number of neighbors of each grid 
### (sum of num must match with the length of adj)

## pre-calculate factorial
# lgamma(x) is equivalent to the natural-log of the factorial = log((x-1)!)
# So, lgamma(x+1) = log(x!)
logFactorial = lgamma((1:gsBreaks[gs_class_n]) + 1)

####################################################################################

path <- 'D:/UserData/KUDrive/@Projects/HKK_Line_Transect/'
source(paste0(path, 'model-detection-function-grsize-landscape-CAR-Kumar2021-22072025.R'))

####################################################################################

constants <- list(I = tran_n, J = dist_class_n,
                  K = gs_class_n, gs_max = gs_max, 
                  L = grid_n, n_covar = ncol(covar),
                  distBreaks = distBreaks, gsBreaks = gsBreaks,
                  adj = adj, num = num, njoin = njoin,
                  logFactorial=logFactorial)
data <- list(y = y_matrix, covar = covar, propM = propM)
inits <- list(muc = runif(1, 1, gs_max), sigma0 = runif(1, 2, 5), p = runif(1, 1, 5), 
              tau = runif(1,0.5,1), beta = rnorm(ncol(covar), 1,2), 
              b_spatial = runif(grid_n, 0, 0.1), w = rep(1, ncol(covar)))

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
#distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data)

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = c("sigma", "p", "muc", "gs_k", 
                                                      "sigma0", "pi", 'lam', 'tr_lam', 
                                                      'beta', 'tau', 'b_spatial', 'w', 
                                                      'AGS', 'ABUND', 'TOTAL_ABUND'))
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

###################################################################################################
# Run the MCMC
t_start <- Sys.time() # start time
samples_chains <- runMCMC(Cmcmc, niter = 1000, nburnin = 500, thin = 2, nchains = 2)
(t_elapse <- Sys.time() - t_start) # time elapse
###summary(samples_chains)

### SAVE MODEL TO RDS files (So you dont have to run the same model again)
saveRDS(samples_chains, 'D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\Model_1.rds')
samples_chains <- readRDS('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\Model_1.rds')

####################################################################################################
################################ Post training ####################################################
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

### merge multiple chain into 1
samples <- do.call(rbind, samples_chains)

### Detection function parameters
sigma_all <- samples[,grep("^sigma", colnames(samples))]
sigma_k <- apply(sigma_all[,1:gs_class_n], 2, median)
print(sigma_k)
p <- median(samples[,'p'])
p

# multinomial pi
pi_val <- apply(samples[,grep("^pi", colnames(samples))], 2, median)
pi <- matrix(pi_val, nrow=gs_class_n, byrow = F)
pi
sum(pi)

### Visualize pi = joint detection prob
plot(pi[1,], type='n', xlab='Distance Class', ylab='Group Size', main='Multinomial Prob')
for(i in 1:nrow(pi)){
  lines(pi[i,], lty=i)
}
legend('topright', legend=1:nrow(pi), title = 'Group size class', lty=1:nrow(pi))

# calculate gx (conditional detection prob given group size class)
x_seq <- seq(0, dist_limit)
hp <- hist(data_sub$P.dist, main='Histogram & Detection function', xlab='Distance')
for(i in 1:length(sigma_k)){
  gx <- exp((-1*x_seq^2)/(2*sigma_k[i]^2))
  lines(x_seq, gx*max(hp$counts), col=i, lwd=2)
}
legend('topright', legend=1:length(sigma_k), title = 'Group size class', col=1:length(sigma_k), lty=1)

# Beta
beta <- samples[,grep("^beta", colnames(samples))]
colnames(beta) <- colnames(covar)
boxplot(beta, outline=F)
abline(h=0, lty=3)

# Indicator
w <- samples[,grep("^w", colnames(samples))]
colnames(w) <- colnames(covar)
w_prop <- colSums(w)/nrow(w)
barplot(w_prop)

# plot beta that w has value at least 50% of 1
beta_w <- beta[,w_prop > 0.5]
boxplot(beta_w, outline=F)
abline(h=0, lty=3)

# Site-specific spatial random effect
b_spatial <- samples[,grep("^b_spatial", colnames(samples))]
# median of spat random effect
bspat_med <- apply(b_spatial, 2, median, na.rm = TRUE)
hist(bspat_med)
# SD of spat random effect
bspat_sd <- apply(b_spatial, 2, sd, na.rm = TRUE)

# Site-specific Lambda (Abundance)
#lam <- samples[,grep("^lam", colnames(samples))]
abund <- samples[,grep("^ABUND", colnames(samples))]
colnames(abund) <- 1:grid_n
# boxplot of abundance of grid 1 to 30
boxplot(abund[,1:30], ylim=c(0,50), outline=F)
# get 5% median and 95% from posterior of abundance
abun_med <- apply(abund, 2, median, na.rm = TRUE)
hist(abun_med, breaks=80)
# SD of  posterior abundance
abun_sd <- apply(abund, 2, sd, na.rm = TRUE)

# Inject back to shp file
library(terra)
# median abundance
poly$abun_med <- abun_med
poly$bspat_med <- bspat_med
plot(poly["abun_med"], main='Median grid-level abundance (lambda) \n Fixed + Random Effects',
     nbreaks = 20, breaks='quantile')
plot(poly["bspat_med"], main='CAR \n Median spatial random effects',
     breaks='quantile')
# SD abundance
poly$abun_sd <- abun_sd
poly$bspat_sd <- bspat_sd
plot(poly["abun_sd"], main='SD grid-level abundance (lambda) \n Fixed + Random Effects',
     nbreaks = 20, breaks='quantile')
plot(poly["bspat_sd"], main='CAR \n SD spatial random effects',
     breaks='quantile')

# Total abundance
total_abund <- samples[,grep("^TOTAL_ABUND", colnames(samples))]
hist(total_abund, breaks=20)
abline(v=median(total_abund), lty=2)
quantile(total_abund, c(0.05, 0.25, 0.5, 0.75, 0.95))
