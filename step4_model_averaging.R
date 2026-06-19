# step4_model_averaging.R
library(nimble)
library(dplyr)
library(coda)
library(sf)
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

cat("Starting Step 4: Model Averaging and Prediction Mapping\n")

if (!file.exists("Results/Final_Model_Comparison.csv")) {
  stop("Could not find Results/Final_Model_Comparison.csv. Please wait for step23 to finish.")
}
final_results <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)

species_list <- unique(final_results$Species)
all_species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK")

# Create output dir for tiffs
if(!dir.exists("Results/Maps")) {
  dir.create("Results/Maps", recursive=TRUE)
}

for (sp_name in species_list) {
  cat("\n========================================================================\n")
  cat("SPECIES: ", sp_name, "\n")
  cat("========================================================================\n")
  
  sp_code <- all_species_codes[[sp_name]]
  
  # Select models to average (Delta_wAIC <= 2)
  sp_res <- final_results %>% filter(Species == sp_name, Delta_wAIC <= 2)
  if (nrow(sp_res) == 0) {
    cat("No models found for", sp_name, "\n")
    next
  }
  
  # Recalculate normalized weights just among the selected models
  sp_res$NormWeight <- exp(-0.5 * sp_res$Delta_wAIC) / sum(exp(-0.5 * sp_res$Delta_wAIC))
  cat(sprintf("Averaging %d models...\n", nrow(sp_res)))
  
  # Setup species data
  species <- sp_code
  nrep <- 10
  dist_limit <- 100
  data_sub_sp <- data_tr %>% filter(Species == sp_code)
  data_tr_sp <- data_tr
  data_tr_sp$P.dist[data_tr_sp$P.dist > dist_limit] <- dist_limit
  gs_max_sp <- max(data_sub_sp$Gz.sz)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  dist_class_n <- 5
  
  # Initialize the model-averaged abundance array
  avg_abund <- rep(0, nrow(poly))
  
  for (i in 1:nrow(sp_res)) {
    model_type <- sp_res$Type[i]
    covars_str <- sp_res$Covariates[i]
    weight <- sp_res$NormWeight[i]
    
    cat("  Running Model", i, "/", nrow(sp_res), "-", model_type, "-", covars_str, "(Weight:", round(weight, 3), ")\n")
    
    # 1. Prepare covariate mask
    w_mask <- rep(0, n_covar)
    if (covars_str != "Intercept-only") {
      active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
      w_mask[covar_names %in% active_vars] <- 1
    }
    
    # 2. Build NIMBLE model
    garbage_out <- capture.output({
      source('@data_prepare_011025.R', local=TRUE)
      if (model_type == "Spatial") {
        source('CAR-Kumar2021-IND.R', local=TRUE)
      } else {
        source('CAR-Kumar2021-NoSpatial.R', local=TRUE)
      }
      
      distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
      distanceModel$w <- w_mask
      distanceModel$calculate()
      
      # We only need ABUND, which contains z[l]*AGS for all 3013 grid cells
      mcmcConf <- configureMCMC(distanceModel, monitors = c("ABUND"))
      
      # We must remove the sampler for the indicator variables 'w' since we fix them!
      mcmcConf$removeSamplers('w')
      
      distanceMCMC <- buildMCMC(mcmcConf)
      CdistanceModel <- compileNimble(distanceModel)
      Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)
      
      # Set `w` mask manually into compiled model
      CdistanceModel$w <- w_mask
      CdistanceModel$calculate()
      
      # Run MCMC
      samples <- runMCMC(Cmcmc, niter = 50000, nburnin = 30000, thin = 2, nchains = 1, setSeed = 999)
    })
    
    # columns are ABUND[1], ABUND[2]... ABUND[3013]
    pred_mean <- colMeans(samples)
    
    # Add to the weighted average
    avg_abund <- avg_abund + (pred_mean * weight)
  }
  
  # 3. Export to GeoTIFF
  cat("  Exporting Map...\n")
  sp_poly <- poly
  sp_poly$Pred_Abund <- avg_abund
  
  vect_poly <- vect(sp_poly)
  # Create a raster template with 1km resolution based on extent of poly
  template <- rast(ext(vect_poly), resolution = c(1000, 1000), crs = crs(vect_poly))
  
  # Rasterize the polygons
  r_abund <- rasterize(vect_poly, template, field = "Pred_Abund")
  
  # Write to tif
  out_name <- paste0("Results/Maps/Abundance_", sp_code, ".tif")
  writeRaster(r_abund, out_name, overwrite=TRUE)
  cat("  Saved:", out_name, "\n")
}

cat("Finished all model averaging and map generation!\n")
