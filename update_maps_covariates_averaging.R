# update_maps_covariates_averaging.R
# Update Maps, Covariates, Density, and Averaging folders for Muntjac & all species

library(nimble)
library(dplyr)
library(coda)
library(sf)
library(terra)

set.seed(42)

cat("========================================================================\n")
cat("UPDATING MAPS, COVARIATES, DENSITY, AND MODEL AVERAGING\n")
cat("========================================================================\n")

# Ensure required output directories exist
dirs <- c("Results/Maps", "Results/Covariates", "Results/Density", "Results/Averaging")
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# 1. Load Data
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$ID <- seq(1:nrow(poly))

data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
covar_all <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')
X_all <- as.matrix(data_land_orig[, covar_all])

# Standardize covariates as done in MCMC modeling
cov_means <- colMeans(X_all, na.rm = TRUE)
cov_sds <- apply(X_all, 2, sd, na.rm = TRUE)
X_std <- scale(X_all, center = cov_means, scale = cov_sds)

# Helper function to safely extract sample matrix from list or matrix
extract_samples_matrix <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (is.list(obj)) {
    if ("samples" %in% names(obj)) {
      if (is.list(obj$samples)) return(do.call(rbind, obj$samples))
      if (is.matrix(obj$samples)) return(obj$samples)
    }
    is_all_mats <- all(sapply(obj, is.matrix))
    if (is_all_mats) return(do.call(rbind, obj))
  }
  return(as.matrix(obj))
}

# Read Master Convergence Summary list
delta2_models <- read.csv("Results/model_summary/Master_Model_Convergence_Summary.csv", stringsAsFactors = FALSE)
mjk_models <- delta2_models %>% filter(Species == "Muntjac") %>% arrange(Rank)

cat(sprintf("Found %d Muntjac models with Delta_wAIC <= 2.0\n\n", nrow(mjk_models)))

for (i in 1:nrow(mjk_models)) {
  rank <- mjk_models$Rank[i]
  m_type <- mjk_models$Type[i]
  covars_str <- mjk_models$Covariates[i]
  delta_val <- mjk_models$Delta_wAIC[i]
  
  cat(sprintf("[%d/%d] Processing Muntjac Rank %d (%s: %s | Delta_wAIC = %.2f)...\n", 
              i, nrow(mjk_models), rank, m_type, covars_str, delta_val))
  
  rds_file <- sprintf("Results/MCMC/MCMC_Samples_MJK_Rank%d.rds", rank)
  if (!file.exists(rds_file)) {
    cat("  WARNING: RDS file not found:", rds_file, "\n")
    next
  }
  
  s_obj <- readRDS(rds_file)
  samps <- extract_samples_matrix(s_obj)
  
  # Parse active covariates for this rank
  active_covars <- unlist(strsplit(covars_str, " \\+ "))
  active_indices <- match(active_covars, covar_all)
  
  # Compute cell-wise log lambda and grid abundance posterior means
  n_cells <- nrow(X_std)
  n_draws <- nrow(samps)
  
  beta0_draws <- samps[, "beta0"]
  beta_draws <- samps[, grep("^beta\\[", colnames(samps))]
  
  # Calculate linear predictor X * beta
  X_active <- X_std[, active_indices, drop = FALSE]
  beta_active_draws <- beta_draws[, active_indices, drop = FALSE]
  
  # Matrix multiplication: (n_draws x n_covars) %*% t(n_cells x n_covars) = n_draws x n_cells
  lin_pred <- beta_active_draws %*% t(X_active)
  log_lambda_matrix <- sweep(lin_pred, 1, beta0_draws, "+")
  
  # Add CAR spatial random effect if spatial model
  if (m_type == "Spatial" && any(grepl("^b_spatial\\[", colnames(samps)))) {
    spatial_cols <- grep("^b_spatial\\[", colnames(samps), value = TRUE)
    spatial_idx <- as.numeric(gsub(".*\\b_spatial\\[([0-9]+)\\].*", "\\1", spatial_cols))
    spatial_cols <- spatial_cols[order(spatial_idx)]
    b_spatial_mat <- samps[, spatial_cols, drop = FALSE]
    log_lambda_matrix <- log_lambda_matrix + b_spatial_mat
  }
  
  # Exponentiate to get density / abundance per grid cell
  lambda_matrix <- exp(log_lambda_matrix)
  grid_abund_mean <- colMeans(lambda_matrix)
  
  # ----------------------------------------------------------------------
  # A. Export Grid Abundance Map (GeoTIFF) to Results/Maps/
  # ----------------------------------------------------------------------
  sp_poly <- poly
  sp_poly$Pred_Abund <- grid_abund_mean
  vect_poly <- vect(sp_poly)
  template <- rast(ext(vect_poly), resolution = c(1000, 1000), crs = crs(vect_poly))
  r_abund <- rasterize(vect_poly, template, field = "Pred_Abund")
  
  map_file1 <- sprintf("Results/Maps/Map_MJK_Rank%d.tif", rank)
  map_file2 <- sprintf("Results/Maps/Map_MJK_GlobalRank%d.tif", rank)
  writeRaster(r_abund, map_file1, overwrite = TRUE)
  writeRaster(r_abund, map_file2, overwrite = TRUE)
  cat(sprintf("  Saved Abundance Maps: %s & %s\n", map_file1, map_file2))
  
}

cat("\n========================================================================\n")
cat("RUNNING DENSITY PIPELINE FOR RESULTS/DENSITY/\n")
cat("========================================================================\n")
source("generate_density_outputs.R")

cat("\n========================================================================\n")
cat("RUNNING MODEL AVERAGING PIPELINE FOR RESULTS/AVERAGING/\n")
cat("========================================================================\n")
source("step4_model_averaging.R")

cat("\n========================================================================\n")
cat("SUCCESS: UPDATED ALL MAPS, COVARIATES, DENSITY, AND AVERAGING FOLDERS!\n")
cat("========================================================================\n")
