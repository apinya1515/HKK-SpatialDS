library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)

# Load data
data_tr <- read.table('line_data.txt', sep='\t', header=T)
poly <- st_read('shp/HKK1sqkmGrid.shp')
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')

# Initialization
species <- "BTG"
nrep <- 10
dist_limit <- 100
data_tr$P.dist[data_tr$P.dist > dist_limit ] <- dist_limit
gsBreaks <- c(1, 2, 3, 4, 8)
dist_class_n <- 5

# Data preparation
source('@data_prepare_011025.R')

# Model file
source('CAR-Kumar2021-IND.R')

cat("Building NIMBLE spatial model...\n")
distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)

# Configure MCMC with minimal tracked variables
tracked_var_opt <- c("beta0", "beta", "w")
mcmcConf <- configureMCMC(distanceModel, monitors = tracked_var_opt, enableWAIC = TRUE)
mcmcConf$removeSamplers('w')

cat("Building MCMC...\n")
distanceMCMC <- buildMCMC(mcmcConf)

cat("Compiling model (C++)...\n")
CdistanceModel <- compileNimble(distanceModel)

cat("Compiling MCMC (C++)...\n")
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

# Run a 1,000 iteration test
cat("Running 1,000 iterations (optimized monitors)...\n")
t_start <- Sys.time()
samples_chains <- runMCMC(Cmcmc, niter = 1000, nburnin = 500, thin = 1, nchains = 1, WAIC=TRUE, setSeed = 999)
t_elapse <- Sys.time() - t_start

cat(sprintf("Time for 1,000 iterations: %.2f seconds\n", as.numeric(t_elapse, units="secs")))
cat(sprintf("Estimated time for 50,000 iterations: %.2f minutes\n", (as.numeric(t_elapse, units="secs") * 50) / 60))
