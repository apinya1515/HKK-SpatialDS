# generate_car_maps.R
# Generates spatial maps for CAR random effects (b_spatial) across all spatial models (Delta_wAIC <= 2)
# Saves outputs to CAR_map/ and Results/CAR_map/

library(sf)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(coda)

cat("========================================================================\n")
cat("GENERATING CAR SPATIAL RANDOM EFFECT MAPS FOR ALL SPATIAL MODELS\n")
cat("========================================================================\n")

car_dir_root <- "CAR_map"
car_dir_res <- "Results/CAR_map"

if (!dir.exists(car_dir_root)) dir.create(car_dir_root, recursive = TRUE)
if (!dir.exists(car_dir_res)) dir.create(car_dir_res, recursive = TRUE)

mcmc_dir <- "Results/MCMC"
post_dir <- "Results/Posteriors"
delta_dir <- "Results/Delta2_Models"

poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$grid_id <- 1:nrow(poly)

species_info <- list(
  BTG = "Banteng",
  SBR = "Sambar deer",
  GAR = "Gaur",
  MJK = "Muntjac",
  PIG = "Wild boar"
)

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

# Load model selection results
if (file.exists("Results/Final_Model_Comparison.csv")) {
  final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  spatial_models <- final_df %>% filter(Delta_wAIC <= 2 & Type == "Spatial") %>% arrange(Species, Rank)
} else if (file.exists("Results/Delta2_Models/Models_Delta2_Summary.csv")) {
  final_df <- read.csv("Results/Delta2_Models/Models_Delta2_Summary.csv", stringsAsFactors = FALSE)
  spatial_models <- final_df %>% filter(Type == "Spatial") %>% arrange(Species, Rank)
} else {
  stop("Model summary CSV not found.")
}

cat(sprintf("Found %d Spatial models with Delta_wAIC <= 2 across all species.\n\n", nrow(spatial_models)))

extract_samples_matrix <- function(s_obj) {
  if (is.matrix(s_obj)) return(s_obj)
  if (is.list(s_obj)) {
    if ("samples" %in% names(s_obj)) {
      if (is.list(s_obj$samples)) return(do.call(rbind, s_obj$samples))
      if (is.matrix(s_obj$samples)) return(s_obj$samples)
    }
    is_all_matrices <- all(sapply(s_obj, is.matrix))
    if (is_all_matrices) return(do.call(rbind, s_obj))
  }
  return(as.matrix(s_obj))
}

rank1_car_df <- list()

for (i in 1:nrow(spatial_models)) {
  sp_name <- spatial_models$Species[i]
  sp_code <- species_codes[[sp_name]]
  rank <- spatial_models$Rank[i]
  covars_str <- spatial_models$Covariates[i]
  delta_val <- spatial_models$Delta_wAIC[i]
  
  cat(sprintf("Processing CAR Map for [%s Rank %d]: %s (Delta_wAIC = %.2f)...\n", sp_code, rank, covars_str, delta_val))
  
  candidates <- c(
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(delta_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)),
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
  )
  
  found_rds <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      found_rds <- cand
      break
    }
  }
  
  if (is.null(found_rds)) {
    cat(sprintf("  WARNING: MCMC file for %s Rank %d not found. Skipping...\n", sp_code, rank))
    next
  }
  
  samples_matrix <- tryCatch({
    s_obj <- readRDS(found_rds)
    extract_samples_matrix(s_obj)
  }, error = function(e) {
    cat(sprintf("  ERROR loading %s: %s\n", found_rds, e$message))
    NULL
  })
  
  if (is.null(samples_matrix)) next
  
  spatial_cols <- grep("^(b_spatial|spatial_z)\\[", colnames(samples_matrix), value = TRUE)
  
  if (length(spatial_cols) == 0) {
    cat(sprintf("  WARNING: No b_spatial columns in %s. Skipping CAR map...\n", found_rds))
    next
  }
  
  spatial_indices <- as.numeric(gsub(".*\\b_spatial\\[([0-9]+)\\].*", "\\1", spatial_cols))
  if (any(is.na(spatial_indices))) {
    spatial_indices <- as.numeric(gsub(".*\\[([0-9]+)\\].*", "\\1", spatial_cols))
  }
  spatial_cols <- spatial_cols[order(spatial_indices)]
  
  b_spatial <- samples_matrix[, spatial_cols, drop = FALSE]
  
  bspat_median <- apply(b_spatial, 2, median)
  bspat_sd     <- apply(b_spatial, 2, sd)
  bspat_rr     <- apply(b_spatial, 2, function(x) median(exp(x)))
  
  poly_sp <- poly
  poly_sp$CAR_Median <- bspat_median[1:nrow(poly_sp)]
  poly_sp$CAR_SD     <- bspat_sd[1:nrow(poly_sp)]
  poly_sp$CAR_RR     <- bspat_rr[1:nrow(poly_sp)]
  
  if (rank == 1 || !sp_code %in% names(rank1_car_df)) {
    rank1_car_df[[sp_code]] <- data.frame(
      grid_id = 1:nrow(poly_sp),
      Species = sp_name,
      Species_Code = sp_code,
      CAR_Median = poly_sp$CAR_Median,
      CAR_SD = poly_sp$CAR_SD,
      CAR_RR = poly_sp$CAR_RR
    )
  }
  
  # Save CSV of CAR summaries per cell
  car_csv <- data.frame(
    grid_id = 1:nrow(poly_sp),
    Species = sp_name,
    Species_Code = sp_code,
    Rank = rank,
    Covariates = covars_str,
    CAR_Median = poly_sp$CAR_Median,
    CAR_SD = poly_sp$CAR_SD,
    CAR_RelativeRisk = poly_sp$CAR_RR
  )
  
  csv_file1 <- file.path(car_dir_root, sprintf("CAR_Summary_%s_Rank%d.csv", sp_code, rank))
  csv_file2 <- file.path(car_dir_res, sprintf("CAR_Summary_%s_Rank%d.csv", sp_code, rank))
  write.csv(car_csv, csv_file1, row.names = FALSE)
  write.csv(car_csv, csv_file2, row.names = FALSE)
  
  # Limit color range for high contrast
  med_lim <- max(abs(quantile(poly_sp$CAR_Median, c(0.01, 0.99), na.rm=TRUE)))
  if (med_lim == 0 || is.na(med_lim)) med_lim <- 1.0
  
  p1 <- ggplot(poly_sp) +
    geom_sf(aes(fill = CAR_Median), color = NA) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "#FFFFBF", high = "#D7191C", 
                         midpoint = 0, limits = c(-med_lim, med_lim), oob = scales::squish,
                         name = "b_spatial\n(Median)") +
    theme_minimal() +
    labs(title = sprintf("CAR Random Effect: %s (Rank %d)", sp_name, rank),
         subtitle = sprintf("Covariates: %s | Delta_wAIC = %.2f", covars_str, delta_val)) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  p2 <- ggplot(poly_sp) +
    geom_sf(aes(fill = CAR_SD), color = NA) +
    scale_fill_viridis_c(option = "magma", name = "CAR SD\n(Uncertainty)") +
    theme_minimal() +
    labs(title = "Uncertainty (Posterior SD)", subtitle = "Standard Deviation of b_spatial") +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  p3 <- ggplot(poly_sp) +
    geom_sf(aes(fill = CAR_RR), color = NA) +
    scale_fill_viridis_c(option = "viridis", name = "Relative Risk\nexp(b_spatial)") +
    theme_minimal() +
    labs(title = "Relative Risk Multiplier", subtitle = "exp(b_spatial) Abundance Multiplier") +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  combined_map <- grid.arrange(p1, p2, p3, ncol = 3)
  
  map_file1 <- file.path(car_dir_root, sprintf("CAR_Map_%s_Rank%d.png", sp_code, rank))
  map_file2 <- file.path(car_dir_res, sprintf("CAR_Map_%s_Rank%d.png", sp_code, rank))
  
  ggsave(map_file1, combined_map, width = 15, height = 6, dpi = 300)
  ggsave(map_file2, combined_map, width = 15, height = 6, dpi = 300)
  
  cat(sprintf("  Saved CAR map: %s\n", map_file1))
}

# Multi-species CAR Comparison Map
if (length(rank1_car_df) > 0) {
  cat("\nGenerating Multi-Species CAR Comparative Map...\n")
  all_sp_car <- do.call(rbind, rank1_car_df)
  poly_multi <- poly %>% left_join(all_sp_car, by = "grid_id")
  
  p_multi <- ggplot(poly_multi) +
    geom_sf(aes(fill = CAR_Median), color = NA) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "#FFFFBF", high = "#D7191C", 
                         midpoint = 0, name = "CAR Effect\n(b_spatial)") +
    facet_wrap(~Species, ncol = 3) +
    theme_minimal() +
    labs(title = "Multi-Species CAR Spatial Random Effects Comparison",
         subtitle = "Huai Kha Khaeng Wildlife Sanctuary (1 km² Grid)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  multi_file1 <- file.path(car_dir_root, "CAR_Map_All_Species_Comparison.png")
  multi_file2 <- file.path(car_dir_res, "CAR_Map_All_Species_Comparison.png")
  
  ggsave(multi_file1, p_multi, width = 14, height = 10, dpi = 300)
  ggsave(multi_file2, p_multi, width = 14, height = 10, dpi = 300)
  
  cat(sprintf("  Saved Multi-Species CAR map: %s\n", multi_file1))
}

cat("\n========================================================================\n")
cat("CAR MAP GENERATION COMPLETE!\n")
cat("========================================================================\n")
