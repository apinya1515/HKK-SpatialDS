# generate_density_outputs.R
# Generate 3-Panel Posterior Histograms: Cluster Density, Individual Density, and Cluster Size (AGS)
# Output Directory: Results/Density/

library(nimble)
library(dplyr)
library(ggplot2)
library(gridExtra)

set.seed(42)

cat("========================================================================\n")
cat("GENERATING 3-PANEL DENSITY & CLUSTER SIZE POSTERIOR HISTOGRAMS\n")
cat("Panels: 1) Cluster Density | 2) Individual Density | 3) Cluster Size (AGS)\n")
cat("========================================================================\n")

density_dir <- "Results/Density"
if (!dir.exists(density_dir)) dir.create(density_dir, recursive = TRUE)

# Remove legacy Results/Posteriors folder if present
if (dir.exists("Results/Posteriors")) {
  cat("Removing legacy Results/Posteriors/ directory...\n")
  unlink("Results/Posteriors", recursive = TRUE)
}

study_area <- 3011 # km2

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

master_df <- read.csv("Results/model_summary/Master_Model_Convergence_Summary.csv", stringsAsFactors = FALSE)

cat(sprintf("Processing 3-panel histograms for %d top models across 5 species...\n\n", nrow(master_df)))

density_summary_list <- list()
all_density_draws_list <- list()

for (i in 1:nrow(master_df)) {
  sp_name <- master_df$Species[i]
  sp_code <- master_df$Species_Code[i]
  rank <- master_df$Rank[i]
  m_type <- master_df$Type[i]
  covars_str <- master_df$Covariates[i]
  delta_val <- master_df$Delta_wAIC[i]
  
  cat(sprintf("[%d/%d] Generating 3-panel Density Histograms for %s Rank %d (%s: %s)...\n", 
              i, nrow(master_df), sp_code, rank, m_type, covars_str))
  
  rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s_Rank%d.rds", sp_code, rank)
  if (!file.exists(rds_file)) {
    rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s.rds", sp_code)
  }
  if (!file.exists(rds_file)) {
    cat("  WARNING: MCMC RDS file not found for", sp_code, "Rank", rank, "\n")
    next
  }
  
  s_obj <- readRDS(rds_file)
  samps <- extract_samples_matrix(s_obj)
  
  tot_abund <- samps[, "TOTAL_ABUND"]
  
  if ("AGS" %in% colnames(samps)) {
    ags_draws <- samps[, "AGS"]
  } else {
    ags_draws <- rep(1, length(tot_abund))
  }
  
  # Calculate Density Metrics
  cluster_counts <- tot_abund / ags_draws
  d_cluster_km2 <- cluster_counts / study_area
  d_ind_km2 <- tot_abund / study_area
  
  d_cluster_100km2 <- d_cluster_km2 * 100
  d_ind_100km2 <- d_ind_km2 * 100
  
  # Summary Statistics
  summary_row <- data.frame(
    Species = sp_name,
    Species_Code = sp_code,
    Rank = rank,
    Type = m_type,
    Covariates = covars_str,
    Delta_wAIC = delta_val,
    Cluster_Density_Mean_km2 = mean(d_cluster_km2),
    Cluster_Density_Median_km2 = median(d_cluster_km2),
    Cluster_Density_SD_km2 = sd(d_cluster_km2),
    Cluster_Density_95_CI_km2 = sprintf("[%.3f, %.3f]", quantile(d_cluster_km2, 0.025), quantile(d_cluster_km2, 0.975)),
    Cluster_Density_Mean_100km2 = mean(d_cluster_100km2),
    Cluster_Density_Median_100km2 = median(d_cluster_100km2),
    Cluster_Density_95_CI_100km2 = sprintf("[%.1f, %.1f]", quantile(d_cluster_100km2, 0.025), quantile(d_cluster_100km2, 0.975)),
    Ind_Density_Mean_km2 = mean(d_ind_km2),
    Ind_Density_Median_km2 = median(d_ind_km2),
    Ind_Density_SD_km2 = sd(d_ind_km2),
    Ind_Density_95_CI_km2 = sprintf("[%.3f, %.3f]", quantile(d_ind_km2, 0.025), quantile(d_ind_km2, 0.975)),
    Ind_Density_Mean_100km2 = mean(d_ind_100km2),
    Ind_Density_Median_100km2 = median(d_ind_100km2),
    Ind_Density_95_CI_100km2 = sprintf("[%.1f, %.1f]", quantile(d_ind_100km2, 0.025), quantile(d_ind_100km2, 0.975)),
    Cluster_Size_Mean = mean(ags_draws),
    Cluster_Size_Median = median(ags_draws),
    Cluster_Size_SD = sd(ags_draws),
    Cluster_Size_95_CI = sprintf("[%.2f, %.2f]", quantile(ags_draws, 0.025), quantile(ags_draws, 0.975)),
    stringsAsFactors = FALSE
  )
  density_summary_list[[length(density_summary_list) + 1]] <- summary_row
  
  # Subsample for comparative plotting
  set.seed(100 + i)
  sub_idx <- sample(1:length(tot_abund), min(1000, length(tot_abund)))
  all_density_draws_list[[length(all_density_draws_list) + 1]] <- data.frame(
    Species = sp_name,
    Species_Code = sp_code,
    Rank = rank,
    Cluster_Density_km2 = d_cluster_km2[sub_idx],
    Ind_Density_km2 = d_ind_km2[sub_idx],
    Cluster_Size = ags_draws[sub_idx],
    stringsAsFactors = FALSE
  )
  
  # Helper to make unsquished histogram with staggered text labels
  make_annotated_hist <- function(values, title_text, xlab_text, fill_color, is_int = FALSE) {
    df_val <- data.frame(Val = values)
    m_val <- mean(values)
    med_val <- median(values)
    q2.5_val <- quantile(values, 0.025)
    q97.5_val <- quantile(values, 0.975)
    q99_val <- quantile(values, 0.99)
    
    # Cap x-axis to prevent squishing from extreme outliers
    x_min <- max(0, q2.5_val * 0.7)
    x_max <- q99_val * 1.25
    if (x_max == x_min) x_max <- x_min + 1.0
    
    df_filtered <- df_val %>% filter(Val <= x_max)
    
    fmt_str <- if (is_int) "%.2f" else "%.3f"
    
    ggplot(df_filtered, aes(x = Val)) +
      geom_histogram(bins = 45, fill = fill_color, color = "white", alpha = 0.85) +
      # Mean Line (Red Solid)
      geom_vline(xintercept = m_val, color = "#D9534F", linewidth = 1.2, linetype = "solid") +
      # Median Line (Blue Dashed)
      geom_vline(xintercept = med_val, color = "#0275D8", linewidth = 1.2, linetype = "dashed") +
      # 95% CI Lines (Dark Red Dotted)
      geom_vline(xintercept = q2.5_val, color = "#A94442", linewidth = 1.0, linetype = "dotted") +
      geom_vline(xintercept = q97.5_val, color = "#A94442", linewidth = 1.0, linetype = "dotted") +
      # Staggered Inward Text Labels
      annotate("text", x = m_val, y = Inf, label = sprintf(paste0("Mean: ", fmt_str), m_val), 
               color = "#D9534F", fontface = "bold", vjust = 2.2, hjust = -0.15, size = 3.3) +
      annotate("text", x = med_val, y = Inf, label = sprintf(paste0("Median: ", fmt_str), med_val), 
               color = "#0275D8", fontface = "bold", vjust = 4.2, hjust = 1.15, size = 3.3) +
      annotate("text", x = q2.5_val, y = Inf, label = sprintf(paste0("2.5%%: ", fmt_str), q2.5_val), 
               color = "#A94442", fontface = "bold", vjust = 6.2, hjust = 1.15, size = 3.1) +
      annotate("text", x = q97.5_val, y = Inf, label = sprintf(paste0("97.5%%: ", fmt_str), q97.5_val), 
               color = "#A94442", fontface = "bold", vjust = 6.2, hjust = 1.15, size = 3.1) +
      scale_x_continuous(limits = c(x_min, x_max)) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11.5),
            plot.subtitle = element_text(size = 8.5, color = "gray30"),
            axis.title = element_text(face = "bold", size = 10),
            panel.grid.minor = element_blank(),
            plot.margin = margin(10, 15, 10, 15)) +
      labs(title = title_text,
           subtitle = sprintf(paste0("Mean: ", fmt_str, " | Med: ", fmt_str, " | 95%% CI: [", fmt_str, ", ", fmt_str, "]"), 
                              m_val, med_val, q2.5_val, q97.5_val),
           x = xlab_text, y = "Frequency")
  }
  
  p1 <- make_annotated_hist(d_cluster_km2, 
                            sprintf("Cluster Density (D_cluster): %s Rank %d", sp_name, rank),
                            "Cluster Density (clusters / km²)", "#48C9B0")
  
  p2 <- make_annotated_hist(d_ind_km2, 
                            sprintf("Individual Density (D_ind): %s Rank %d", sp_name, rank),
                            "Individual Density (individuals / km²)", "#5DADE2")
  
  p3 <- make_annotated_hist(ags_draws, 
                            sprintf("Cluster Size (AGS): %s Rank %d", sp_name, rank),
                            "Cluster Size (mean individuals / cluster)", "#F5B041", is_int = TRUE)
  
  p_comb <- grid.arrange(p1, p2, p3, ncol = 3)
  
  out_png <- sprintf("Results/Density/Density_Hist_%s_Rank%d.png", sp_code, rank)
  ggsave(out_png, p_comb, width = 16, height = 5.5, dpi = 300)
  cat(sprintf("  Saved 3-Panel Density & Cluster Size Histogram: %s\n", out_png))
}

# Save Summary Table
master_density_df <- do.call(rbind, density_summary_list)
write.csv(master_density_df, "Results/Density/Species_Density_Cluster_and_Individual_Summary.csv", row.names = FALSE)
cat("\n  Saved Density Summary CSV: Results/Density/Species_Density_Cluster_and_Individual_Summary.csv\n")

# ------------------------------------------------------------------------
# Master Multi-Species Comparison Figure for Density & Cluster Size
# ------------------------------------------------------------------------
cat("\nGenerating Master Multi-Species Density & Cluster Size Comparison Figure...\n")

if (length(all_density_draws_list) > 0) {
  draws_df <- do.call(rbind, all_density_draws_list)
  rank1_draws <- draws_df %>% filter(Rank == 1)
  
  p_master_clust <- ggplot(rank1_draws, aes(x = Species, y = Cluster_Density_km2, fill = Species)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.5) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(face = "bold", angle = 20, hjust = 1),
          legend.position = "none") +
    labs(title = "Cluster Density (D_cluster)", x = "Species", y = "Clusters / km²")
  
  p_master_ind <- ggplot(rank1_draws, aes(x = Species, y = Ind_Density_km2, fill = Species)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.5) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(face = "bold", angle = 20, hjust = 1),
          legend.position = "none") +
    labs(title = "Individual Density (D_individual)", x = "Species", y = "Individuals / km²")
  
  p_master_ags <- ggplot(rank1_draws, aes(x = Species, y = Cluster_Size, fill = Species)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, color = "black", linewidth = 0.5) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(face = "bold", angle = 20, hjust = 1),
          legend.position = "none") +
    labs(title = "Cluster Size (AGS)", x = "Species", y = "Individuals / Cluster")
  
  p_master_comp <- grid.arrange(p_master_clust, p_master_ind, p_master_ags, ncol = 3)
  master_density_png <- "Results/Density/Density_Hist_All_Species_Comparison.png"
  ggsave(master_density_png, p_master_comp, width = 16, height = 5.5, dpi = 300)
  cat(sprintf("  Saved Master 3-Panel Density & Cluster Size Comparison Figure: %s\n", master_density_png))
}

cat("\n========================================================================\n")
cat("SUCCESS: 3-PANEL DENSITY & CLUSTER SIZE POSTERIOR HISTOGRAMS COMPLETED!\n")
cat("========================================================================\n")
