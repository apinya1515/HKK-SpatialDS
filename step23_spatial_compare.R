# step23_spatial_compare.R
# 1. Reads the NoSpatial model selection results
# 2. Extracts top 5 models
# 3. Runs the CAR-Kumar2021-IND.R (spatial) model for those covariate combinations
# 4. Computes a global Delta WAIC and Weight across both NoSpatial and Spatial models
# 5. Computes quantiles of TOTAL_ABUND
# 6. Generates raster map of grid 50% abundance
#
# ARCHITECTURE: Sequential single-process execution.
#   - Compile NIMBLE model ONCE per species, then loop through covariate masks.
#   - This avoids the massive RAM overhead of parallel compilation (~8 GB for 3 workers)
#     and the std::bad_alloc crashes that result from it.

library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)
library(terra)

options(nimbleVerbose = FALSE)

# 1. Load Data
data_tr <- read.table('line_data.txt', sep='\t', header=T)
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')

covar_names <- colnames(data_land)[-1]
n_covar <- length(covar_names)

cat("Starting Step 2 & 3: Spatial Model Evaluation & Final Comparison\n")

if(!dir.exists("Results/Maps")) dir.create("Results/Maps", recursive=TRUE)
if(!dir.exists("Results/Posteriors")) dir.create("Results/Posteriors", recursive=TRUE)

# Read NoSpatial results
if (!file.exists("Results/Model_Selection_Summary.csv")) {
  stop("Could not find Results/Model_Selection_Summary.csv. Please run model_selection.R first.")
}
nospatial_results <- read.csv("Results/Model_Selection_Summary.csv", stringsAsFactors = FALSE)

# We will collect the final combined results here
final_all_species_results <- list()

species_list <- unique(nospatial_results$Species)
all_species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK")

# Check for already-completed species (crash recovery)
completed_species <- c()
for (sp_check in names(all_species_codes)) {
  sp_code_check <- all_species_codes[[sp_check]]
  intermediate_file <- sprintf("Results/Intermediate_Comparison_%s.csv", sp_code_check)
  if (file.exists(intermediate_file)) {
    cat(sprintf("Found intermediate results for %s, loading from file.\n", sp_check))
    final_all_species_results[[sp_code_check]] <- read.csv(intermediate_file, stringsAsFactors = FALSE)
    completed_species <- c(completed_species, sp_check)
  }
}

for (sp_name in species_list) {
  # Skip already-completed species
  if (sp_name %in% completed_species) {
    cat(sprintf("\nSkipping %s (already completed from previous run).\n", sp_name))
    next
  }

  cat("\n========================================================================\n")
  cat("SPECIES: ", sp_name, "\n")
  cat("========================================================================\n")
  
  sp_code <- all_species_codes[[sp_name]]
  
  # Filter NoSpatial results for this species
  sp_nospatial <- nospatial_results %>% filter(Species == sp_name) %>% arrange(Rank)
  sp_nospatial$Type <- "NoSpatial"
  
  # Select top 5 models to run spatial formulation
  target_models <- head(sp_nospatial, 5)
  
  cat(sprintf("Evaluating Spatial variants for the top %d NoSpatial models.\n", nrow(target_models)))
  
  if(nrow(target_models) == 0) {
    # If none, just use NoSpatial for final comparison
    final_all_species_results[[sp_code]] <- sp_nospatial
    next
  }
  
  # Determine w_masks for these models
  w_masks_to_run <- list()
  for (i in 1:nrow(target_models)) {
    covars_str <- target_models$Covariates[i]
    w <- rep(0, n_covar)
    if (covars_str != "Intercept-only") {
      active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
      w[covar_names %in% active_vars] <- 1
    }
    w_masks_to_run[[i]] <- list(mask = w, name = covars_str)
  }
  
  # Configure species parameters
  species <- sp_code
  nrep <- 10
  dist_limit <- 100
  data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
  data_sub_sp <- data_tr %>% filter(Species == sp_code)
  gs_max_sp <- max(data_sub_sp$Gz.sz)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  dist_class_n <- 5
  
  # ---- LOOP THROUGH COVARIATE MASKS (RECOMPILE EACH TO AVOID MEMORY LEAK) ----
  spatial_results_list <- list()
  
  for (i in seq_along(w_masks_to_run)) {
    w_item <- w_masks_to_run[[i]]
    rank <- i
    covars_str <- w_item$name
    
    # ---- MODEL-LEVEL CRASH RECOVERY ----
    # Check if this specific model already has saved RDS results
    rds_file <- sprintf("Results/Posteriors/Samples_%s_Rank%d.rds", sp_code, rank)
    if (file.exists(rds_file)) {
      cat(sprintf("[%s] Model %d/%d: %s -> SKIPPING (found %s)\n", format(Sys.time(), "%H:%M:%S"), i, length(w_masks_to_run), covars_str, rds_file))
      flush.console()
      
      # Reload saved results to extract summary stats
      tryCatch({
        saved_chains <- readRDS(rds_file)
        samps <- saved_chains$samples
        total_abund_samples <- samps[, "TOTAL_ABUND"]
        abund_cols <- grep("^ABUND\\[", colnames(samps))
        abund_med <- apply(samps[, abund_cols], 2, median)
        quants <- quantile(total_abund_samples, probs = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975))
        
        df_row <- data.frame(
          Species = sp_name, Rank = NA, Covariates = covars_str,
          wAIC = saved_chains$WAIC$WAIC, Delta_wAIC = NA, Weight = NA,
          lppd = saved_chains$WAIC$lppd, pWAIC = saved_chains$WAIC$pWAIC,
          Type = "Spatial",
          CI_2.5 = quants["2.5%"], CI_5 = quants["5%"], CI_25 = quants["25%"],
          CI_50 = quants["50%"], CI_75 = quants["75%"], CI_95 = quants["95%"],
          CI_97.5 = quants["97.5%"], stringsAsFactors = FALSE
        )
        colnames(df_row)[10:16] <- c("CI_2.5", "CI_5", "CI_25", "CI_50", "CI_75", "CI_95", "CI_97.5")
        spatial_results_list[[i]] <- df_row
        
        cat(sprintf("[%s]   -> Loaded: WAIC = %.2f, TOTAL_ABUND median = %.1f\n",
                    format(Sys.time(), "%H:%M:%S"), saved_chains$WAIC$WAIC, median(total_abund_samples)))
        rm(saved_chains, samps, total_abund_samples, abund_med, quants)
        gc()
      }, error = function(e) {
        cat(sprintf("[%s]   -> WARNING: Could not reload %s: %s. Will re-run.\n", format(Sys.time(), "%H:%M:%S"), rds_file, e$message))
        file.remove(rds_file)  # Remove corrupt file so it gets re-run
      })
      
      # If we successfully loaded, skip to next model
      if (!is.null(spatial_results_list[[i]])) next
    }
    
    cat(sprintf("[%s] Model %d/%d: %s\n", format(Sys.time(), "%H:%M:%S"), i, length(w_masks_to_run), covars_str))
    cat(sprintf("[%s]   Compiling NIMBLE model fresh (avoids memory leak)...\n", format(Sys.time(), "%H:%M:%S")))
    flush.console()
    
    res <- tryCatch({
      # ---- FRESH COMPILATION FOR EACH MODEL ----
      garbage_out <- capture.output({
        source('@data_prepare_011025.R')
        source('CAR-Kumar2021-IND.R')
        
        distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
        mcmcConf <- configureMCMC(distanceModel, monitors = c("beta0", "beta", "w", "ABUND", "TOTAL_ABUND"), enableWAIC = TRUE)
        mcmcConf$removeSamplers('w')
        distanceMCMC <- buildMCMC(mcmcConf)
        CdistanceModel <- compileNimble(distanceModel)
        Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)
      })
      rm(garbage_out)
      
      cat(sprintf("[%s]   Compilation done. Starting MCMC...\n", format(Sys.time(), "%H:%M:%S")))
      flush.console()
      
      # Set the covariate mask and recalculate
      CdistanceModel$w <- w_item$mask
      CdistanceModel$calculate()
      
      garbage_mcmc <- capture.output({
        samples_chains <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000, thin = 2, nchains = 1, WAIC = TRUE,
                                  setSeed = as.integer(Sys.time()) + i)
      })
      rm(garbage_mcmc)
      
      samps <- samples_chains$samples
      
      total_abund_samples <- samps[, "TOTAL_ABUND"]
      abund_cols <- grep("^ABUND\\[", colnames(samps))
      abund_med <- apply(samps[, abund_cols], 2, median)
      
      result <- list(
        Covariates = covars_str,
        WAIC = samples_chains$WAIC$WAIC,
        lppd = samples_chains$WAIC$lppd,
        pWAIC = samples_chains$WAIC$pWAIC,
        total_abund_samples = total_abund_samples,
        abund_med = abund_med,
        samples_chains = samples_chains,  # Keep for saveRDS
        error = NULL
      )
      
      # ---- CLEAN UP NIMBLE OBJECTS IMMEDIATELY ----
      rm(CdistanceModel, Cmcmc, distanceModel, distanceMCMC, mcmcConf, samps)
      nimble::clearCompiled(distanceModelCode)
      gc()
      
      result
    }, error = function(e) {
      # Clean up on error too
      tryCatch({
        rm(CdistanceModel, Cmcmc, distanceModel, distanceMCMC, mcmcConf, envir = parent.frame())
      }, error = function(e2) {})
      gc()
      
      list(
        Covariates = covars_str,
        WAIC = Inf,
        lppd = NA,
        pWAIC = NA,
        total_abund_samples = NULL,
        abund_med = NULL,
        samples_chains = NULL,
        error = e$message
      )
    })
    
    if (is.null(res$error)) {
      cat(sprintf("[%s]   -> WAIC = %.2f, TOTAL_ABUND median = %.1f\n", 
                  format(Sys.time(), "%H:%M:%S"), res$WAIC, median(res$total_abund_samples)))
      
      # Calculate CI
      quants <- quantile(res$total_abund_samples, probs = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975))
      
      # Save Posterior plot
      jpeg(sprintf("Results/Posteriors/Posterior_%s_Rank%d.jpg", sp_code, rank), width=800, height=600)
      hist(res$total_abund_samples, breaks=50, main=sprintf("Total Abundance Posterior: %s (Rank %d)", sp_name, rank), xlab="Total Abundance", col="lightgray", border="white")
      
      # Median
      abline(v=quants["50%"], col="blue", lwd=2)
      text(quants["50%"], par("usr")[4]*0.9, paste("Median:", round(quants["50%"], 1)), col="blue", pos=4)
      
      # 90% CI
      abline(v=quants["5%"], col="green", lwd=2, lty=2)
      abline(v=quants["95%"], col="green", lwd=2, lty=2)
      text(quants["5%"], par("usr")[4]*0.8, paste("90% L:", round(quants["5%"], 1)), col="darkgreen", pos=2)
      text(quants["95%"], par("usr")[4]*0.8, paste("90% U:", round(quants["95%"], 1)), col="darkgreen", pos=4)
      
      # 95% CI
      abline(v=quants["2.5%"], col="red", lwd=2, lty=2)
      abline(v=quants["97.5%"], col="red", lwd=2, lty=2)
      text(quants["2.5%"], par("usr")[4]*0.7, paste("95% L:", round(quants["2.5%"], 1)), col="red", pos=2)
      text(quants["97.5%"], par("usr")[4]*0.7, paste("95% U:", round(quants["97.5%"], 1)), col="red", pos=4)
      
      dev.off()
      
      # Save Map
      poly_map <- poly
      poly_map$abund_med <- res$abund_med
      vect_poly <- vect(poly_map)
      template <- rast(ext(vect_poly), resolution = c(1000, 1000), crs = crs(vect_poly))
      r_abund <- rasterize(vect_poly, template, field = "abund_med")
      writeRaster(r_abund, sprintf("Results/Maps/Map_%s_Rank%d.tif", sp_code, rank), overwrite=TRUE)
      
      # Save full MCMC posterior samples to RDS for later analysis
      saveRDS(res$samples_chains, sprintf("Results/Posteriors/Samples_%s_Rank%d.rds", sp_code, rank))
      
      df_row <- data.frame(
        Species = sp_name,
        Rank = NA,
        Covariates = res$Covariates,
        wAIC = res$WAIC,
        Delta_wAIC = NA,
        Weight = NA,
        lppd = res$lppd,
        pWAIC = res$pWAIC,
        Type = "Spatial",
        CI_2.5 = quants["2.5%"],
        CI_5 = quants["5%"],
        CI_25 = quants["25%"],
        CI_50 = quants["50%"],
        CI_75 = quants["75%"],
        CI_95 = quants["95%"],
        CI_97.5 = quants["97.5%"],
        stringsAsFactors = FALSE
      )
      colnames(df_row)[10:16] <- c("CI_2.5", "CI_5", "CI_25", "CI_50", "CI_75", "CI_95", "CI_97.5")
      spatial_results_list[[i]] <- df_row
    } else {
      cat(sprintf("[%s]   -> FAILED: %s\n", format(Sys.time(), "%H:%M:%S"), res$error))
    }
    
    # ---- AGGRESSIVE CLEANUP BETWEEN MODELS ----
    rm(res)
    gc()
    cat(sprintf("[%s]   Memory after cleanup: %.0f MB used\n", format(Sys.time(), "%H:%M:%S"), sum(gc()[,2])))
    flush.console()
  }

  
  # ---- COMBINE NOSPATIAL + SPATIAL RESULTS ----
  if (length(spatial_results_list) > 0) {
    spatial_df <- do.call(rbind, spatial_results_list)
    
    # For NoSpatial, add empty CI columns so we can rbind
    sp_nospatial$CI_2.5 <- NA; sp_nospatial$CI_5 <- NA; sp_nospatial$CI_25 <- NA
    sp_nospatial$CI_50 <- NA; sp_nospatial$CI_75 <- NA; sp_nospatial$CI_95 <- NA; sp_nospatial$CI_97.5 <- NA
    
    # Ensure column order matches
    common_cols <- c("Species", "Rank", "Covariates", "wAIC", "Delta_wAIC", "Weight", "lppd", "pWAIC", "Type", "CI_2.5", "CI_5", "CI_25", "CI_50", "CI_75", "CI_95", "CI_97.5")
    sp_nospatial <- sp_nospatial[, common_cols]
    spatial_df <- spatial_df[, common_cols]
    
    # Combine NoSpatial and Spatial
    combined_df <- rbind(sp_nospatial, spatial_df)
    
    # Recalculate global Rank, Delta wAIC, and Weight
    combined_df <- combined_df %>% arrange(wAIC)
    best_waic <- combined_df$wAIC[1]
    combined_df$Delta_wAIC <- combined_df$wAIC - best_waic
    combined_df$Weight <- exp(-0.5 * combined_df$Delta_wAIC) / sum(exp(-0.5 * combined_df$Delta_wAIC))
    combined_df$Rank <- 1:nrow(combined_df)
    
    final_all_species_results[[sp_code]] <- combined_df
    
    cat(sprintf("\nGlobal Top 5 models for %s:\n", sp_name))
    print(head(combined_df, 5))
    
    # Save intermediate progress (enables crash recovery)
    write.csv(combined_df, sprintf("Results/Intermediate_Comparison_%s.csv", sp_code), row.names = FALSE)
    cat(sprintf("[%s] Saved intermediate results for %s.\n", format(Sys.time(), "%H:%M:%S"), sp_name))
  } else {
    cat(sprintf("\nAll spatial models failed for %s. Using NoSpatial.\n", sp_name))
    final_all_species_results[[sp_code]] <- sp_nospatial
  }
  
  # Force garbage collection to free RAM before next species
  # (NIMBLE objects are cleaned up per-model now)
  gc()
}

cat("\n========================================================================\n")
cat("FINAL COMPARISON COMPLETE\n")
cat("========================================================================\n")

final_summary_combined <- do.call(rbind, final_all_species_results)
write.csv(final_summary_combined, "Results/Final_Model_Comparison.csv", row.names = FALSE)
cat("Saved to 'Results/Final_Model_Comparison.csv'.\n")
