# step4_model_averaging.R
library(dplyr)
library(coda)
library(sf)
library(terra)

options(nimbleVerbose = FALSE)

# 1. Load Data
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$ID <- seq(1:nrow(poly))

cat("Starting Step 4: Model Averaging via MCMC Sample Pooling\n")

# Try to find Results/Final_Model_Comparison.csv or Intermediate files
results_files <- list.files("Results", pattern = "Comparison", full.names = TRUE)
if (length(results_files) == 0) {
  stop("No model comparison results found. Please wait for step23 to produce output.")
}

# Load all results and combine
final_results_list <- lapply(results_files, function(f) {
  if (grepl("csv$", f)) read.csv(f, stringsAsFactors = FALSE) else NULL
})
final_results <- bind_rows(final_results_list) %>% distinct()

species_list <- unique(final_results$Species)
all_species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK")

# Create output dir for tiffs
if(!dir.exists("Results/Maps")) {
  dir.create("Results/Maps", recursive=TRUE)
}

# Number of total pooled samples we want for final calculation
TOTAL_SAMPLES <- 10000

for (sp_name in species_list) {
  cat("\n========================================================================\n")
  cat("SPECIES: ", sp_name, "\n")
  cat("========================================================================\n")
  
  sp_code <- all_species_codes[[sp_name]]
  
  # Select models to average (Delta_wAIC <= 2)
  sp_res <- final_results %>% filter(Species == sp_name, Delta_wAIC <= 2) %>% arrange(Rank)
  
  # We can only average Spatial models because NoSpatial models do not have ABUND stored
  sp_res <- sp_res %>% filter(Type == "Spatial")
  
  if (nrow(sp_res) == 0) {
    cat("No Spatial models with Delta_wAIC <= 2 found for", sp_name, "\n")
    next
  }
  
  # Recalculate normalized weights
  sp_res$NormWeight <- exp(-0.5 * sp_res$Delta_wAIC) / sum(exp(-0.5 * sp_res$Delta_wAIC))
  cat(sprintf("Averaging %d Spatial models...\n", nrow(sp_res)))
  
  # Initialize lists for pooled samples
  pooled_abund_list <- list()
  pooled_total_list <- list()
  
  for (i in 1:nrow(sp_res)) {
    model_rank <- sp_res$Rank[i]
    covars_str <- sp_res$Covariates[i]
    weight <- sp_res$NormWeight[i]
    
    rds_file <- sprintf("Results/Posteriors/Samples_%s_Rank%d.rds", sp_code, model_rank)
    if (!file.exists(rds_file)) {
      cat(sprintf("  WARNING: Could not find MCMC samples for Model Rank %d. Skipping...\n", model_rank))
      next
    }
    
    cat("  Loading Model Rank", model_rank, "-", covars_str, "(Weight:", round(weight, 3), ")\n")
    
    s_obj <- readRDS(rds_file)
    samps <- s_obj$samples
    
    # Calculate how many samples to draw from this model
    n_draw <- round(TOTAL_SAMPLES * weight)
    
    if (n_draw > 0) {
      # Randomly sample rows from MCMC
      draw_idx <- sample(1:nrow(samps), n_draw, replace = TRUE)
      
      # Extract ABUND columns
      abund_cols <- grep("^ABUND\\[", colnames(samps))
      drawn_abund <- samps[draw_idx, abund_cols]
      pooled_abund_list[[i]] <- drawn_abund
      
      # Extract TOTAL_ABUND
      drawn_total <- samps[draw_idx, "TOTAL_ABUND"]
      pooled_total_list[[i]] <- drawn_total
    }
  }
  
  # Combine pooled samples
  pooled_abund_matrix <- do.call(rbind, pooled_abund_list)
  pooled_total_vector <- unlist(pooled_total_list)
  
  actual_samples <- nrow(pooled_abund_matrix)
  cat(sprintf("\n  Pooled a total of %d MCMC samples across models.\n", actual_samples))
  
  if (actual_samples == 0) {
    cat("  No samples pooled. Skipping species.\n")
    next
  }
  
  # Calculate grid-wise Mean and SD
  cat("  Calculating Grid-wise Variance...\n")
  grid_mean <- apply(pooled_abund_matrix, 2, mean)
  grid_sd   <- apply(pooled_abund_matrix, 2, sd)
  
  # Calculate Total Abundance stats
  total_mean <- mean(pooled_total_vector)
  total_median <- median(pooled_total_vector)
  total_ci <- quantile(pooled_total_vector, probs = c(0.025, 0.975))
  
  cat("  --------------------------------------\n")
  cat(sprintf("  TOTAL ABUNDANCE BMA ESTIMATES:\n"))
  cat(sprintf("  Mean:   %.1f\n", total_mean))
  cat(sprintf("  Median: %.1f\n", total_median))
  cat(sprintf("  95%% CI: %.1f - %.1f\n", total_ci[1], total_ci[2]))
  cat("  --------------------------------------\n")
  
  # Export to GeoTIFF
  cat("  Exporting Maps...\n")
  sp_poly <- poly
  sp_poly$Pred_Abund <- grid_mean
  sp_poly$Pred_SD <- grid_sd
  
  vect_poly <- vect(sp_poly)
  template <- rast(ext(vect_poly), resolution = c(1000, 1000), crs = crs(vect_poly))
  
  r_abund <- rasterize(vect_poly, template, field = "Pred_Abund")
  out_name_abund <- paste0("Results/Maps/Abundance_", sp_code, ".tif")
  writeRaster(r_abund, out_name_abund, overwrite=TRUE)
  cat("  Saved:", out_name_abund, "\n")
  
  r_sd <- rasterize(vect_poly, template, field = "Pred_SD")
  out_name_sd <- paste0("Results/Maps/SD_Abundance_", sp_code, ".tif")
  writeRaster(r_sd, out_name_sd, overwrite=TRUE)
  cat("  Saved:", out_name_sd, "\n")
  
  # Free up memory
  rm(pooled_abund_matrix, pooled_abund_list, samps, s_obj)
  gc()
}

cat("Finished all model averaging and map generation!\n")
