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
all_species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

# Rename existing _Rank files to _GlobalRank files based on NoSpatial covariates
cat("\nChecking for legacy _Rank files to rename to _GlobalRank...\n")
for (sp_name in species_list) {
  sp_code <- all_species_codes[[sp_name]]
  sp_nospatial <- final_results %>% filter(Species == sp_name, Type == "NoSpatial") %>% arrange(wAIC)
  sp_spatial <- final_results %>% filter(Species == sp_name, Type == "Spatial")
  
  if (nrow(sp_nospatial) == 0 || nrow(sp_spatial) == 0) next
  
  for (i in 1:nrow(sp_spatial)) {
    global_rank <- sp_spatial$Rank[i]
    covars_str <- sp_spatial$Covariates[i]
    
    nospatial_idx <- which(sp_nospatial$Covariates == covars_str)
    if (length(nospatial_idx) == 0) next
    nospatial_rank <- nospatial_idx[1]
    
    # Define old and new filenames
    old_rds <- sprintf("Results/Posteriors/Samples_%s_Rank%d.rds", sp_code, nospatial_rank)
    new_rds <- sprintf("Results/Posteriors/Samples_%s_GlobalRank%d.rds", sp_code, global_rank)
    if (file.exists(old_rds) && !file.exists(new_rds)) file.rename(old_rds, new_rds)
    
    old_jpg <- sprintf("Results/Posteriors/Posterior_%s_Rank%d.jpg", sp_code, nospatial_rank)
    new_jpg <- sprintf("Results/Posteriors/Posterior_%s_GlobalRank%d.jpg", sp_code, global_rank)
    if (file.exists(old_jpg) && !file.exists(new_jpg)) file.rename(old_jpg, new_jpg)
    
    old_tif <- sprintf("Results/Maps/Map_%s_Rank%d.tif", sp_code, nospatial_rank)
    new_tif <- sprintf("Results/Maps/Map_%s_GlobalRank%d.tif", sp_code, global_rank)
    if (file.exists(old_tif) && !file.exists(new_tif)) file.rename(old_tif, new_tif)
  }
}
cat("Renaming check complete.\n")

# Create output dir for averaging results
if(!dir.exists("Results/Averaging")) {
  dir.create("Results/Averaging", recursive=TRUE)
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
  
  if (nrow(sp_res) == 0) {
    cat("No models with Delta_wAIC <= 2 found for", sp_name, "\n")
    next
  }
  
  has_spatial <- any(sp_res$Type == "Spatial")
  if (!has_spatial) {
    cat("  Note: Only NoSpatial models found. Will skip mapping but will average TOTAL_ABUND.\n")
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
    
    rds_candidates <- c(
      sprintf("Results/MCMC/MCMC_Samples_%s_Rank%d.rds", sp_code, model_rank),
      sprintf("Results/Delta2_Models/Samples_%s_Rank%d.rds", sp_code, model_rank),
      sprintf("Results/Posteriors/Samples_%s_GlobalRank%d.rds", sp_code, model_rank),
      sprintf("Results/Posteriors/Samples_%s_Rank%d.rds", sp_code, model_rank),
      sprintf("Results/MCMC/MCMC_Samples_%s.rds", sp_code)
    )
    rds_file <- NULL
    for (cand in rds_candidates) {
      if (file.exists(cand)) {
        rds_file <- cand
        break
      }
    }
    if (is.null(rds_file)) {
      cat(sprintf("  WARNING: Could not find MCMC samples for %s Global Rank %d. Skipping...\n", sp_code, model_rank))
      next
    }
    
    cat("  Loading Global Rank", model_rank, "-", covars_str, "from", rds_file, "(Weight:", round(weight, 3), ")\n")
    
    s_obj <- readRDS(rds_file)
    samps <- s_obj$samples
    
    # Calculate how many samples to draw from this model
    n_draw <- round(TOTAL_SAMPLES * weight)
    
    if (n_draw > 0) {
      # Randomly sample rows from MCMC
      draw_idx <- sample(1:nrow(samps), n_draw, replace = TRUE)
      
      # Extract ABUND columns if Spatial
      if (sp_res$Type[i] == "Spatial") {
        abund_cols <- grep("^ABUND\\[", colnames(samps))
        drawn_abund <- samps[draw_idx, abund_cols]
        pooled_abund_list[[length(pooled_abund_list) + 1]] <- drawn_abund
      }
      
      # Extract TOTAL_ABUND
      drawn_total <- samps[draw_idx, "TOTAL_ABUND"]
      pooled_total_list[[length(pooled_total_list) + 1]] <- drawn_total
    }
  }
  
  # Combine pooled samples
  pooled_total_vector <- unlist(pooled_total_list)
  
  actual_samples <- length(pooled_total_vector)
  cat(sprintf("\n  Pooled a total of %d MCMC samples across models.\n", actual_samples))
  
  if (actual_samples == 0) {
    cat("  No samples pooled. Skipping species.\n")
    next
  }
  
  if (length(pooled_abund_list) > 0) {
    pooled_abund_matrix <- do.call(rbind, pooled_abund_list)
    cat("  Calculating Grid-wise Variance...\n")
    grid_mean <- apply(pooled_abund_matrix, 2, mean)
    grid_sd   <- apply(pooled_abund_matrix, 2, sd)
    grid_cv   <- ifelse(grid_mean == 0, 0, grid_sd / grid_mean)
  } else {
    grid_mean <- NULL
  }
  
  # Calculate Total Abundance stats
  total_mean <- mean(pooled_total_vector)
  total_median <- median(pooled_total_vector)
  total_ci <- quantile(pooled_total_vector, probs = c(0.025, 0.975))
  total_ci_90 <- quantile(pooled_total_vector, probs = c(0.05, 0.95))
  
  cat("  --------------------------------------\n")
  cat(sprintf("  TOTAL ABUNDANCE BMA ESTIMATES:\n"))
  cat(sprintf("  Mean:   %.1f\n", total_mean))
  cat(sprintf("  Median: %.1f\n", total_median))
  cat(sprintf("  90%% CI: %.1f - %.1f\n", total_ci_90[1], total_ci_90[2]))
  cat(sprintf("  95%% CI: %.1f - %.1f\n", total_ci[1], total_ci[2]))
  cat("  --------------------------------------\n")
  
  # Plot Histogram
  cat("  Saving Histogram...\n")
  jpeg(paste0("Results/Averaging/TotalAbundance_Hist_", sp_code, ".jpg"), width=800, height=600)
  x_max <- quantile(pooled_total_vector, probs = 0.975) * 1.15
  plot_data <- pooled_total_vector[pooled_total_vector <= x_max]
  hist(plot_data, breaks=50, main=paste("Model Averaged Total Abundance:", sp_name), 
       xlab="Total Abundance (Truncated at 97.5th %ile * 1.15)", col="lightgray", border="white", xlim=c(min(plot_data), x_max))
  
  # Median line
  abline(v=total_median, col="blue", lwd=2)
  text(total_median, par("usr")[4]*0.9, paste("Median:", round(total_median, 1)), col="blue", pos=4)
  
  # 90% CI
  abline(v=total_ci_90[1], col="green", lwd=2, lty=2)
  abline(v=total_ci_90[2], col="green", lwd=2, lty=2)
  text(total_ci_90[1], par("usr")[4]*0.8, paste("90% L:", round(total_ci_90[1], 1)), col="darkgreen", pos=2)
  text(total_ci_90[2], par("usr")[4]*0.8, paste("90% U:", round(total_ci_90[2], 1)), col="darkgreen", pos=4)
  
  # 95% CI
  abline(v=total_ci[1], col="red", lwd=2, lty=2)
  abline(v=total_ci[2], col="red", lwd=2, lty=2)
  text(total_ci[1], par("usr")[4]*0.7, paste("95% L:", round(total_ci[1], 1)), col="red", pos=2)
  text(total_ci[2], par("usr")[4]*0.7, paste("95% U:", round(total_ci[2], 1)), col="red", pos=4)
  
  dev.off()
  
  if (!is.null(grid_mean)) {
    # Export to GeoTIFF
    cat("  Exporting Maps...\n")
    sp_poly <- poly
    sp_poly$Pred_Abund <- grid_mean
    sp_poly$Pred_SD <- grid_sd
    sp_poly$Pred_CV <- grid_cv
    
    vect_poly <- vect(sp_poly)
    template <- rast(ext(vect_poly), resolution = c(1000, 1000), crs = crs(vect_poly))
    
    r_abund <- rasterize(vect_poly, template, field = "Pred_Abund")
    out_name_abund <- paste0("Results/Averaging/Abundance_", sp_code, ".tif")
    writeRaster(r_abund, out_name_abund, overwrite=TRUE)
    cat("  Saved:", out_name_abund, "\n")
    
    r_sd <- rasterize(vect_poly, template, field = "Pred_SD")
    out_name_sd <- paste0("Results/Averaging/SD_Abundance_", sp_code, ".tif")
    writeRaster(r_sd, out_name_sd, overwrite=TRUE)
    cat("  Saved:", out_name_sd, "\n")
    
    r_cv <- rasterize(vect_poly, template, field = "Pred_CV")
    out_name_cv <- paste0("Results/Averaging/CV_Abundance_", sp_code, ".tif")
    writeRaster(r_cv, out_name_cv, overwrite=TRUE)
    cat("  Saved:", out_name_cv, "\n")
    
    rm(pooled_abund_matrix)
  } else {
    cat("  Skipping map export (No Spatial models in top subset).\n")
  }
  
  # Free up memory
  rm(pooled_abund_list, samps, s_obj)
  gc()
}

cat("Finished all model averaging and map generation!\n")
