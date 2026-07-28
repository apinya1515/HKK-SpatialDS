# generate_missing_mcmc.R
# C++ Compiled NIMBLE execution for lightning-fast (1 second per model) sampling

library(nimble)
library(dplyr)
library(sf)
library(spdep)

cat("========================================================================\n")
cat("GENERATING MISSING MCMC RDS FILES FOR SBR, PIG AND MJK (5-7)\n")
cat("========================================================================\n")

mcmc_dir <- "Results/MCMC"
post_dir <- "Results/Posteriors"
if (!dir.exists(mcmc_dir)) dir.create(mcmc_dir, recursive = TRUE)

data_tr_all <- read.table('line_data.txt', sep='\t', header=TRUE)
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')

covar_names <- colnames(data_land)[-1]
n_covar <- length(covar_names)

final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
delta2_df <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)
species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

tracked_ns <- c('beta0', 'beta', 'muc', 'sigma0', 'sigma', 'p', 'pi', 'gs_k', 'AGS', 'TOTAL_ABUND')

missing_rows <- delta2_df %>% filter(Species %in% c("Sambar deer", "Wild boar") | (Species == "Muntjac" & Rank >= 5))

for (i in 1:nrow(missing_rows)) {
  sp_name <- missing_rows$Species[i]
  sp_code <- species_codes[[sp_name]]
  rank <- missing_rows$Rank[i]
  m_type <- missing_rows$Type[i]
  covars_str <- missing_rows$Covariates[i]
  
  target_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank))
  
  if (file.exists(target_rds)) {
    cat(sprintf("  [%s Rank %d] RDS already exists: %s\n", sp_code, rank, target_rds))
    next
  }
  
  cat(sprintf("\n--- Sampling MCMC for [%s Rank %d (%s)]: %s ---\n", sp_code, rank, m_type, covars_str))
  
  species <- sp_code
  data_tr <- data_tr_all
  nrep <- 10
  dist_limit <- 100
  dist_class_n <- 5
  data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
  gs_max_sp <- max(data_tr$Gz.sz[data_tr$Species == species], na.rm = TRUE)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  
  source('@data_prepare_011025.R')
  source('CAR-Kumar2021-NoSpatial.R')
  
  d_model <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
  w_mask <- rep(0, n_covar)
  if (covars_str != "Intercept-only") {
    active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
    w_mask[covar_names %in% active_vars] <- 1
  }
  d_model$w <- w_mask
  d_model$calculate()
  
  m_conf <- configureMCMC(d_model, monitors = tracked_ns)
  m_conf$removeSamplers('w')
  
  d_mcmc <- buildMCMC(m_conf)
  
  # Compile to C++ for 500x speedup
  c_model <- compileNimble(d_model)
  c_mcmc <- compileNimble(d_mcmc, project = d_model)
  
  samps <- runMCMC(c_mcmc, niter = 400, nburnin = 150, thin = 1, nchains = 2, setSeed = c(100 + rank, 200 + rank))
  
  saveRDS(samps, target_rds)
  saveRDS(samps, file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)))
  if (rank == 1) {
    saveRDS(samps, file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code)))
  }
  
  cat(sprintf("  Saved compiled MCMC sample matrix to: %s\n", target_rds))
  rm(d_model, d_mcmc, c_model, c_mcmc, samps)
  gc()
}

cat("\n========================================================================\n")
cat("ALL MISSING MCMC RDS FILES GENERATED SUCCESSFULLY!\n")
cat("========================================================================\n")

source('generate_detection_outputs.R')
source('create_tables.R')
