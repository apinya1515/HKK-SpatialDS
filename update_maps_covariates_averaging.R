# update_maps_covariates_averaging.R
# Update Maps (all Delta_wAIC <= 2 models), Density, and Averaging folders.
# Maps: 4-band GeoTIFF of INDIVIDUAL abundance per cell (mean, median, SD) + extrapolation flag;
# Results/Maps/Extrapolation_Summary.csv gives the share of abundance in extrapolated cells.

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
X_all[X_all[, "ndvi_cv"] > 0.5, "ndvi_cv"] <- 0.5 # same cap as @data_prepare_011025.R

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

# Read Master Convergence Summary list (all species; previously Muntjac only)
delta2_models <- read.csv("Results/model_summary/Master_Model_Convergence_Summary.csv", stringsAsFactors = FALSE)
map_models <- delta2_models %>% arrange(Species, Rank)

cat(sprintf("Found %d models with Delta_wAIC <= 2.0\n\n", nrow(map_models)))

water_mask <- ifelse(data_land_orig$WA > 0.35, 0, 1) # as in @data_prepare_011025.R
# Covariate range (standardised) over the grid cells crossed by transects
tr_cells <- unique(read.csv("TRidentity.csv")$grid_id)
tr_range <- apply(X_std[tr_cells, , drop = FALSE], 2, range)
written <- character(0)
extrap_rows <- list()

for (i in 1:nrow(map_models)) {
  sp_code <- map_models$Species_Code[i]
  rank <- map_models$Rank[i]
  m_type <- map_models$Type[i]
  covars_str <- map_models$Covariates[i]
  delta_val <- map_models$Delta_wAIC[i]

  cat(sprintf("[%d/%d] Processing %s Rank %d (%s: %s | Delta_wAIC = %.2f)...\n",
              i, nrow(map_models), sp_code, rank, m_type, covars_str, delta_val))

  rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s_Rank%d.rds", sp_code, rank)
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
    # (the old pattern ".*\\b_spatial\\[..." never matched -> NA -> file order was used)
    spatial_idx <- as.numeric(gsub("^.*\\[([0-9]+)\\]$", "\\1", spatial_cols))
    stopifnot(!anyNA(spatial_idx))
    spatial_cols <- spatial_cols[order(spatial_idx)]
    b_spatial_mat <- samps[, spatial_cols, drop = FALSE]
    log_lambda_matrix <- log_lambda_matrix + b_spatial_mat
  }
  
  # Individual abundance per grid cell and draw: groups (exp(lin)) x water mask x AGS.
  # Cell means sum to the posterior mean of TOTAL_ABUND.
  abund_matrix <- sweep(exp(log_lambda_matrix), 2, water_mask, "*") * samps[, "AGS"]

  # Extrapolation: any active covariate outside the range covered by transect cells
  extrap <- as.integer(apply(X_std[, active_indices, drop = FALSE], 1, function(x)
    any(x < tr_range[1, active_indices] | x > tr_range[2, active_indices])))

  # ----------------------------------------------------------------------
  # A. Export 4-band abundance map (GeoTIFF) to Results/Maps/
  #    band 1 abund_mean, 2 abund_median, 3 abund_sd (individuals per 1-km2 cell),
  #    band 4 extrapolated (1 = covariates outside the range sampled by transects)
  # ----------------------------------------------------------------------
  sp_poly <- poly
  sp_poly$abund_mean <- colMeans(abund_matrix)
  sp_poly$abund_median <- apply(abund_matrix, 2, median)
  sp_poly$abund_sd <- apply(abund_matrix, 2, sd)
  sp_poly$extrapolated <- extrap
  vect_poly <- vect(sp_poly)
  template <- rast(ext(vect_poly), resolution = c(1000, 1000), crs = crs(vect_poly))
  bands <- c("abund_mean", "abund_median", "abund_sd", "extrapolated")
  r_abund <- rast(lapply(bands, function(b) rasterize(vect_poly, template, field = b)))
  names(r_abund) <- bands

  map_file2 <- sprintf("Results/Maps/Map_%s_GlobalRank%d.tif", sp_code, rank)
  writeRaster(r_abund, map_file2, overwrite = TRUE, names = bands, datatype = "FLT4S") # float: an integer layer would truncate all bands
  written <- c(written, basename(map_file2))
  if (sp_code == "MJK") { # Muntjac maps have always also been saved under the _Rank name
    map_file1 <- sprintf("Results/Maps/Map_MJK_Rank%d.tif", rank)
    writeRaster(r_abund, map_file1, overwrite = TRUE, names = bands, datatype = "FLT4S") # float: an integer layer would truncate all bands
    written <- c(written, basename(map_file1))
  }
  tot <- rowSums(abund_matrix)
  extrap_rows[[length(extrap_rows) + 1]] <- data.frame(
    Species = map_models$Species[i], Species_Code = sp_code, Rank = rank, Type = m_type, Covariates = covars_str,
    Cells_Extrapolated = sum(extrap), Cells_Total = length(extrap),
    Share_of_Mean_Abundance_in_Extrapolated_Cells = round(sum(sp_poly$abund_mean[extrap == 1]) / mean(tot), 4),
    Total_Abundance_Mean = round(mean(tot), 1), Total_Abundance_Median = round(median(tot), 1),
    Total_Abundance_Median_Excluding_Extrapolated = round(median(rowSums(abund_matrix[, extrap == 0, drop = FALSE])), 1),
    stringsAsFactors = FALSE)
  cat(sprintf("  Saved Abundance Map: %s (%d extrapolated cells)\n", map_file2, sum(extrap)))
  rm(s_obj, samps, lin_pred, log_lambda_matrix, abund_matrix); gc()
}

write.csv(do.call(rbind, extrap_rows), "Results/Maps/Extrapolation_Summary.csv", row.names = FALSE)
# Remove maps of models no longer in the Delta_wAIC < 2 set (e.g. old GlobalRank maps of excluded models)
stale <- setdiff(list.files("Results/Maps", pattern = "^Map_.*\\.tif$"), written)
if (length(stale) > 0) {
  unlink(file.path("Results/Maps", stale))
  cat("Removed maps of models outside the Delta_wAIC < 2 set:", paste(stale, collapse = ", "), "\n")
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
