# run_fast_posteriors.R
# Ultra-fast posterior estimation for all 5 species (BTG, SBR, GAR, MJK, PIG)
# Monitors summary parameters for 1-second MCMC sampling per model!

library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)

cat("========================================================================\n")
cat("ULTRA-FAST COMPLETE MCMC POSTERIOR COMPUTATION FOR ALL SPECIES\n")
cat("========================================================================\n")

mcmc_dir <- "Results/MCMC"
post_dir <- "Results/Posteriors"
delta_dir <- "Results/Delta2_Models"

if (!dir.exists(mcmc_dir)) dir.create(mcmc_dir, recursive = TRUE)
if (!dir.exists(post_dir)) dir.create(post_dir, recursive = TRUE)
if (!dir.exists(delta_dir)) dir.create(delta_dir, recursive = TRUE)

data_tr_all <- read.table('line_data.txt', sep='\t', header=TRUE)
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')

covar_names <- colnames(data_land)[-1]
n_covar <- length(covar_names)

species_info <- list(
  BTG = "Banteng",
  SBR = "Sambar deer",
  GAR = "Gaur",
  MJK = "Muntjac",
  PIG = "Wild boar"
)

if (file.exists("Results/Delta2_Models/Models_Delta2_Summary.csv")) {
  delta2_df <- read.csv("Results/Delta2_Models/Models_Delta2_Summary.csv", stringsAsFactors = FALSE)
} else if (file.exists("Results/Final_Model_Comparison.csv")) {
  final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  delta2_df <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)
}

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

tracked_ns <- c('beta0', 'beta', 'muc', 'sigma0', 'sigma', 'p', 'pi', 'gs_k', 'AGS', 'TOTAL_ABUND')
tracked_sp <- c('beta0', 'beta', 'muc', 'sigma0', 'sigma', 'p', 'pi', 'gs_k', 'AGS', 'sigma_spatial', 'TOTAL_ABUND')

for (sp_code in names(species_info)) {
  sp_name <- species_info[[sp_code]]
  cat("\n========================================================================\n")
  cat(sprintf("PROCESSING SPECIES: %s (%s)\n", sp_name, sp_code))
  cat("========================================================================\n")
  
  species <- sp_code
  data_tr <- data_tr_all
  nrep <- 10
  dist_limit <- 100
  dist_class_n <- 5
  data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
  gs_max_sp <- max(data_tr$Gz.sz[data_tr$Species == species], na.rm = TRUE)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  
  source('@data_prepare_011025.R')
  
  sp_models <- delta2_df %>% filter(Species == sp_name)
  if (nrow(sp_models) == 0) {
    sp_models <- data.frame(Species = sp_name, Rank = 1, Type = "NoSpatial", Covariates = "slope", stringsAsFactors = FALSE)
  }
  
  # Run NoSpatial Model for this species
  cat(sprintf("  Running NoSpatial model for %s...\n", sp_name))
  source('CAR-Kumar2021-NoSpatial.R')
  
  distanceModel_ns <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
  mcmcConf_ns <- configureMCMC(distanceModel_ns, monitors = tracked_ns, enableWAIC = TRUE)
  mcmcConf_ns$removeSamplers('w')
  distanceMCMC_ns <- buildMCMC(mcmcConf_ns)
  
  samples_ns <- runMCMC(distanceMCMC_ns, niter = 1000, nburnin = 300, thin = 1, nchains = 2, setSeed = c(101, 202))
  saveRDS(samples_ns, file.path(mcmc_dir, sprintf("MCMC_Samples_%s_NoSpatial.rds", sp_code)))
  
  rm(distanceModel_ns, distanceMCMC_ns)
  gc()
  
  for (m_idx in 1:nrow(sp_models)) {
    rank <- sp_models$Rank[m_idx]
    m_type <- sp_models$Type[m_idx]
    covars_str <- sp_models$Covariates[m_idx]
    
    target_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank))
    target_sp_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
    
    if (m_type == "NoSpatial") {
      file.copy(file.path(mcmc_dir, sprintf("MCMC_Samples_%s_NoSpatial.rds", sp_code)), target_rds, overwrite = TRUE)
      file.copy(file.path(mcmc_dir, sprintf("MCMC_Samples_%s_NoSpatial.rds", sp_code)), file.path(delta_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)), overwrite = TRUE)
      file.copy(file.path(mcmc_dir, sprintf("MCMC_Samples_%s_NoSpatial.rds", sp_code)), file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)), overwrite = TRUE)
      if (rank == 1) file.copy(file.path(mcmc_dir, sprintf("MCMC_Samples_%s_NoSpatial.rds", sp_code)), target_sp_rds, overwrite = TRUE)
      cat(sprintf("  [Rank %d (%s)] Saved NoSpatial posteriors: %s\n", rank, m_type, target_rds))
      next
    }
    
    cat(sprintf("  [Rank %d (%s)] Sampling Spatial model for: %s\n", rank, m_type, covars_str))
    source('CAR-Kumar2021-IND.R')
    
    distanceModel_sp <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
    w_mask <- rep(0, n_covar)
    if (covars_str != "Intercept-only") {
      active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
      w_mask[covar_names %in% active_vars] <- 1
    }
    distanceModel_sp$w <- w_mask
    distanceModel_sp$calculate()
    
    mcmcConf_sp <- configureMCMC(distanceModel_sp, monitors = tracked_sp, enableWAIC = TRUE)
    mcmcConf_sp$removeSamplers('w')
    distanceMCMC_sp <- buildMCMC(mcmcConf_sp)
    
    samples_sp <- runMCMC(distanceMCMC_sp, niter = 1000, nburnin = 300, thin = 1, nchains = 2, setSeed = c(123 + rank, 456 + rank))
    
    saveRDS(samples_sp, target_rds)
    saveRDS(samples_sp, file.path(delta_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)))
    saveRDS(samples_sp, file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)))
    if (rank == 1) saveRDS(samples_sp, target_sp_rds)
    
    cat(sprintf("    Saved Spatial MCMC samples to: %s\n", target_rds))
    rm(distanceModel_sp, distanceMCMC_sp)
    gc()
  }
}

cat("\n========================================================================\n")
cat("ALL POSTERIORS SAVED! REGENERATING TABLES & DETECTION PLOTS...\n")
cat("========================================================================\n")

source('create_tables.R')
source('generate_detection_outputs.R')
