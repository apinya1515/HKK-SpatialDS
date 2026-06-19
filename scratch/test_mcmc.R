# test_mcmc.R - Minimal test to diagnose MCMC crash
cat("=== Starting MCMC crash test ===\n")

tryCatch({
  library(nimble)
  library(dplyr)
  library(coda)
  library(sf)
  library(spdep)
  library(terra)
  options(nimbleVerbose = FALSE)
  cat("Libraries loaded.\n")

  data_tr <- read.table('line_data.txt', sep='\t', header=T)
  poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
  poly$ID <- seq(1:nrow(poly))
  data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
  data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
  data_prop <- read.csv('TRidentity.csv')
  cat("Data loaded.\n")

  species <- 'BTG'
  nrep <- 10
  dist_limit <- 100
  data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
  gsBreaks <- c(1, 2, 3, 4, 8)
  dist_class_n <- 5

  source('@data_prepare_011025.R')
  cat("Data prepared.\n")
  
  source('CAR-Kumar2021-IND.R')
  cat("Model code loaded.\n")

  cat("Building model...\n"); flush.console()
  distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
  
  mcmcConf <- configureMCMC(distanceModel, monitors = c('beta0', 'beta', 'w', 'ABUND', 'TOTAL_ABUND'), enableWAIC = TRUE)
  mcmcConf$removeSamplers('w')
  distanceMCMC <- buildMCMC(mcmcConf)
  cat("Model built.\n"); flush.console()
  
  cat("Compiling model...\n"); flush.console()
  CdistanceModel <- compileNimble(distanceModel)
  cat("Model compiled.\n"); flush.console()
  
  Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)
  cat("MCMC compiled.\n"); flush.console()

  cat("Setting w mask...\n"); flush.console()
  CdistanceModel$w <- rep(1, ncol(covar))
  CdistanceModel$calculate()
  cat("Calculated.\n"); flush.console()

  cat("Running SHORT MCMC (100 iter)...\n"); flush.console()
  samples <- runMCMC(Cmcmc, niter = 100, nburnin = 10, thin = 1, nchains = 1, WAIC = TRUE, setSeed = 999)
  cat("SUCCESS! Short MCMC completed.\n")
  cat(paste("WAIC:", samples$WAIC$WAIC, "\n"))
  
}, error = function(e) {
  cat(paste("ERROR:", e$message, "\n"))
  traceback()
})
