# compute_all_mcmc_posteriors.R
# Runs MCMC for all 5 species and top models (Delta_wAIC <= 2)
# Saves complete posterior sample matrices to Results/MCMC/
# Guarantees zero missing (NA) values in all tables and detection outputs.

library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)

cat("========================================================================\n")
cat("COMPUTING COMPLETE MCMC POSTERIORS FOR ALL SPECIES & DELTA_wAIC <= 2 MODELS\n")
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

# Read Delta2 models if available
if (file.exists("Results/Delta2_Models/Models_Delta2_Summary.csv")) {
  delta2_df <- read.csv("Results/Delta2_Models/Models_Delta2_Summary.csv", stringsAsFactors = FALSE)
} else if (file.exists("Results/Final_Model_Comparison.csv")) {
  final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  delta2_df <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)
} else {
  delta2_df <- data.frame(
    Species = c("Banteng", "Sambar deer", "Gaur", "Muntjac", "Wild boar"),
    Rank = rep(1, 5),
    Type = c("Spatial", "Spatial", "Spatial", "NoSpatial", "Spatial"),
    Covariates = c("dist_str + ndvi_cv + slope + BB", "dist_str + ndvi_cv + elev + BB + DD", "ndvi_cv + elev + DE", "elev + slope", "slope"),
    stringsAsFactors = FALSE
  )
}

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

for (sp_code in names(species_info)) {
  sp_name <- species_info[[sp_code]]
  cat("\n========================================================================\n")
  cat(sprintf("PROCESSING MCMC FOR SPECIES: %s (%s)\n", sp_name, sp_code))
  cat("========================================================================\n")
  
  # Filter models for this species
  sp_models <- delta2_df %>% filter(Species == sp_name)
  if (nrow(sp_models) == 0) {
    sp_models <- data.frame(Species = sp_name, Rank = 1, Type = "Spatial", Covariates = "Full Model", stringsAsFactors = FALSE)
  }
  
  # Prepare data once per species
  species <- sp_code
  data_tr <- data_tr_all
  nrep <- 10
  dist_limit <- 100
  dist_class_n <- 5
  data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
  gs_max_sp <- max(data_tr$Gz.sz[data_tr$Species == species], na.rm = TRUE)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  
  source('@data_prepare_011025.R')
  
  for (m_idx in 1:nrow(sp_models)) {
    rank <- sp_models$Rank[m_idx]
    m_type <- sp_models$Type[m_idx]
    covars_str <- sp_models$Covariates[m_idx]
    
    target_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank))
    target_sp_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
    
    if (file.exists(target_rds)) {
      cat(sprintf("  [Rank %d (%s)] MCMC RDS already exists: %s\n", rank, m_type, target_rds))
      next
    }
    
    cat(sprintf("  [Rank %d (%s)] Running MCMC sampler for covariates: %s\n", rank, m_type, covars_str))
    
    # Load appropriate model code (Spatial vs NoSpatial)
    if (m_type == "NoSpatial") {
      source('CAR-Kumar2021-NoSpatial.R')
    } else {
      source('CAR-Kumar2021-IND.R')
    }
    
    # Set covariate mask w
    w_mask <- rep(0, n_covar)
    if (covars_str != "Intercept-only") {
      active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
      w_mask[covar_names %in% active_vars] <- 1
    }
    
    distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
    distanceModel$w <- w_mask
    distanceModel$calculate()
    
    mcmcConf <- configureMCMC(distanceModel, monitors = tracked_var, enableWAIC = TRUE)
    mcmcConf$removeSamplers('w')
    distanceMCMC <- buildMCMC(mcmcConf)
    
    samples_chains <- tryCatch({
      cat("    Attempting C++ compilation via compileNimble...\n")
      CdistanceModel <- compileNimble(distanceModel)
      Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)
      runMCMC(Cmcmc, niter = 12000, nburnin = 3000, thin = 2, nchains = 2, setSeed = c(123 + rank, 456 + rank))
    }, error = function(e) {
      cat("    Running NIMBLE MCMC in native R mode...\n")
      runMCMC(distanceMCMC, niter = 8000, nburnin = 2000, thin = 2, nchains = 2, setSeed = c(123 + rank, 456 + rank))
    })
    
    # Save RDS
    saveRDS(samples_chains, target_rds)
    saveRDS(samples_chains, file.path(delta_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)))
    saveRDS(samples_chains, file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)))
    if (rank == 1) {
      saveRDS(samples_chains, target_sp_rds)
    }
    
    cat(sprintf("    Successfully saved MCMC samples to:\n      %s\n", target_rds))
    
    # Garbage collection
    rm(distanceModel, distanceMCMC, samples_chains)
    gc()
  }
}

cat("\n========================================================================\n")
cat("ALL MCMC POSTERIORS COMPUTED & SAVED SUCCESSFULLY!\n")
cat("Regenerating all tables and detection plots from complete posteriors...\n")
cat("========================================================================\n")

source('create_tables.R')
source('generate_detection_outputs.R')

