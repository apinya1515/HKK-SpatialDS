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

###### 30/09/2025
### See model code for changes in model
### Split data preparation path into @data_prepare_[date].R

##### 5/10/2025
### Add replication in calculation of transect-level abundance
### lam_tr = nrep * (prop * lam), so, log(lam_tr) = log(nrep) + log(prop * lam)

library(nimble)
library(dplyr)
library(coda) # for post-training mcmc check 
library(MCMCvis)
library(sf) # for import spatial data
library(spdep) # for defining spatial neighbors (used for CAR)
library(terra) # for raster export and visualization

##############################################################################################
# Transect data (table 1)
data_tr <- read.table('D:\\UserData\\Dropbox\\@KU_pop_course\\day1_exercise\\line_data.txt', 
                      sep='\t', header=T)
#data_tr$Tr.no <- paste0(data_tr$Tr.no,'-',data_tr$Walk.no)

# shapefile of study area grid
poly <- st_read('shp\\HKK1sqkmGrid.shp')
poly$ID <- seq(1:nrow(poly))  # define ID directly to data instead
plot(poly)

# landscape data
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_med', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE', 'HE')]

# transect X landscape data
data_prop <- read.csv('TRidentity.csv')

###########################################################################################
########################## INITIALIZATION ####################################
### Input parameters
species <- 'BTG'
# hist(data_tr[data_tr$Species==species,]$P.dist, main='Perp. Distance', breaks=10)
# hist(data_tr[data_tr$Species==species,]$Gz.sz, main='Grp size')
nrep <- 10 # number of replications for each transect
dist_limit <- 100 # max observed distance
data_tr$P.dist[data_tr$P.dist > dist_limit ] <- dist_limit  # set maximum distance to the limit
gsBreaks <- c(1, 4, 6, 8) # upper breaks group size of each gs classes - c(1, 3) -> 2 classes 1st class = 1, 2nd class >= 2, the maximum group size is 3
dist_class_n <- 5 # number of dist classes
##############################################################################

# Set sources path
path <- 'D:/UserData/KUDrive/@Projects/HKK_Line_Transect/'

### Call data preparation
source('@data_prepare_011025.R')

### *****Call model file*****
## Model with indicator
#source(paste0(path, 'model-detection-function-grsize-landscape-CAR-Kumar2021-30092025.R'))
## Model without indicator
source('CAR-Kumar2021-NoSpatial.R')

####################################################################################

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
# see initialization
distanceModel$initializeInfo()

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = tracked_var, enableWAIC = TRUE)
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

###################################################################################################
# Run the MCMC
t_start <- Sys.time() # start time
samples_chains <- runMCMC(Cmcmc, niter = 30000, nburnin = 3000, thin = 2, nchains = 3, WAIC=TRUE,
                          setSeed = c(999, 111, 555))
(t_elapse <- Sys.time() - t_start) # time elapse
###summary(samples_chains)

##############################################
# #### Continuous running
# Cmcmc$run(niter = 5000) # Run for 5000 iter first
# # Then run another more 20000 iterations
# Cmcmc$run(niter = 20000, reset = FALSE) 
# # Then run another 20000
# Cmcmc$run(niter = 20000, reset = FALSE) 

# as.list(Cmcmc$mvSamples)

################################################

### SAVE MODEL TO RDS files (So you dont have to run the same model again)
saveRDS(samples_chains, 'Model_sample.rds')
samples_chains <- readRDS('Model_sample.rds')
#samples_chains <- readRDS("D:/Model_SBR_HKKcut20000burn15000.rds")
####################################################################################################



####################################################################################################
################################ Post training ####################################################
# Convergence check - use package coda
# Diagnostics
# Convert samples to coda format
mcmc.list <- as.mcmc.list(lapply(samples_chains$samples, as.mcmc))
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

### merge multiple chains into 1
samples <- do.call(rbind, samples_chains$samples)
#samples <- samples_chains[[1]]

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

### Indicator (only for model with indicators)
w <- samples[,grep("^w", colnames(samples))]
colnames(w) <- colnames(covar)
w_prop <- colSums(w)/nrow(w)
barplot(w_prop)
# plot beta that w has value at least 10% of 1
beta_w <- beta[,w_prop > 0.3]
boxplot(beta_w, outline=F)
abline(h=0, lty=3)

## Grid-specific spatial random effect of grid abundance (bspat)
b_spatial <- samples[,grep("^b_spatial", colnames(samples))]
# median of spat random effect
bspat_quant <- apply(b_spatial, 2, function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))
hist(bspat_quant[3,], main='Median of Spatial random effects (CAR)')
hist(bspat_quant[5,], main='Upper of Spatial random effects (CAR)')
# SD of spat random effect
bspat_sd <- apply(b_spatial, 2, sd, na.rm = TRUE)
hist(bspat_sd, main='SD of Spatial random effects (CAR)')
# visualization of spatial random effect
poly$bspat_med <- bspat_quant[3,] # Median random spatial effects
plot(poly["bspat_med"], main='Median of spatial random effects (CAR)', breaks='quantile')

### Grid-specific fixed effect of grid group-abundance (z)
z <- samples[,grep("^z", colnames(samples))]
# median of fixed effect
z_quant <- apply(z, 2, function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))
hist(z_quant[3,], main='Median of grid-specific abundance', breaks=100)
hist(z_quant[2,], main='Lower of grid-specific abundance', breaks=100)
hist(z_quant[4,], main='Upper of grid-specific abundance', breaks=100)
hist(z_quant, main='z', breaks=200, xlim=c(0,50))
# SD of fixed effect
z_sd <- apply(z, 2, sd, na.rm = TRUE)
hist(z_sd, main='SD of grid-specific abundance', breaks=50)
# raster visualization for z
poly$z_med <- z_quant[3,] # median of grid abundance
poly$z_upper <- z_quant[4,] # upper of grid abundance
poly$z_lower <- z_quant[2,] # lower of grid abundance
poly$z_sd <- z_sd # SD of fixed effect of grid abundance
plot(poly["z_med"], main='Median of grid-group-abundnace', breaks='quantile')
plot(poly["z_lower"], main='Lower of grid-group-abundnace', breaks='quantile')
plot(poly["z_upper"], main='Upper of grid-group-abundnace', breaks='quantile')
plot(poly["z_sd"], main='SD of grid-group-abundnace', breaks='quantile')

# Transect-specific addtional abundance (alpha)
alpha <- samples[,grep("^alpha", colnames(samples))]
hist(alpha, breaks=50, main='Site-specific additional abundance (alpha)')

# Transect-specific proportionaed grid-abundance (Z)
Z <- samples[,grep("^Z", colnames(samples))]
hist(Z, breaks=100, main='ransect-specific proportionaed grid-abundance (Z)')

### Transect-specific abundance (lam)
lam <- samples[,grep("^lam", colnames(samples))]
boxplot(lam)

### Average group size (AGS)
avg_gs <- samples[,grep("^AGS", colnames(samples))]
hist(avg_gs, main='Median of Average group size (AGS)', breaks=50)

### Grid-specific individual Abundance (ABUND = z * AGS)
abund <- samples[,grep("^ABUND", colnames(samples))]
colnames(abund) <- 1:grid_n
# get 5% median and 95% from posterior of abundance
abund_med <- apply(abund, 2, median, na.rm = TRUE)
abund_upper <- apply(abund, 2, function(x) quantile(x, 0.95))
abund_lower <- apply(abund, 2, function(x) quantile(x, 0.05))
hist(abund_med, breaks=100, xlim=c(0, 100), probability=T, main='Indiv abundance Histogram')
hist(abund_upper, breaks=100, add=T, col=rgb(1,0,0,0.5), probability=T)
hist(abund_lower, breaks=100, add=T, col=rgb(0,1,0,0.5), probability=T)
# SD of  posterior abundance
abund_sd <- apply(abund, 2, sd, na.rm = TRUE)
hist(abund_sd, main='SD of grid-level individual abundance', breaks=100)
# raster visualization for abundance
poly$abund_med <- abund_med # median abundance
poly$abund_upper <- abund_upper
poly$abund_lower <- abund_lower
plot(poly["abund_med"], main='Median grid-level indiv-abundance')
plot(poly["abund_upper"], main='Upper(95%) grid-level indiv-abundance')
plot(poly["abund_lower"], main='Lower(5%) grid-level indiv-abundance')
# SD abundance
poly$abund_sd <- abund_sd
poly$abund_lowsd <- abund_sd < 50
plot(poly["abund_sd"], main='SD grid-level indiv-abundance', nbreaks = 100, breaks='quantile')
plot(poly["abund_lowsd"], main='Abundance - Low SD grid (SD < 100)')
# # subset prediction only low sd
# abund_med_lowsd <- abund[,poly$abund_lowsd]
# dim(abund_lowsd)
# abund_lowsd_quant <- apply(abund_lowsd, 2, function(x) quantile(x, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)))
# abund_lowsd_quant_sum <- apply(abund_lowsd_quant, 1, sum)
# abund_lowsd_quant_sum 

# Median and mean density (abundance per grid) across landscape (not across chain)
abund_med_land <- apply(abund, 1, median, na.rm = TRUE)
hist(abund_med_land, breaks=100)
density_quant <- quantile(abund_med_land, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
density_quant
abline(v=density_quant[1], lty=2);abline(v=density_quant[7], lty=2)
abline(v=density_quant[4])

# Total abundance (from within-model derived-quantities)
total_abund <- samples[,grep("^TOTAL_ABUND", colnames(samples))]
hist(total_abund[total_abund < 100000], breaks=200, main='Total Abundance')
abline(v=median(total_abund), lty=2)
quantile(total_abund, c(0.05, 0.25, 0.5, 0.75, 0.95))
abline(v=quantile(total_abund, c(0.05, 0.95)), col='red', lty=2)

# Calculate total abundace post hoc from grid-wise abundance
sum(abun_med)
sum(abun_upper)
sum(abun_lower)
