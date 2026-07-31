# generate_covariate_boxplots.R
# Generate Covariate Beta Coefficient Boxplots (Covariates Only) with Jitter, Red Mean Marks, 2.5% & 97.5% Quantiles

library(nimble)
library(dplyr)
library(ggplot2)

set.seed(42)

cat("========================================================================\n")
cat("GENERATING COVARIATE BETA BOXPLOTS (COVARIATES ONLY) WITH MARKS & QUANTILES\n")
cat("========================================================================\n")

if (!dir.exists("Results/Covariates")) {
  dir.create("Results/Covariates", recursive = TRUE)
}

covar_names <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')
species_map <- c("BTG" = "Banteng", "GAR" = "Gaur", "MJK" = "Muntjac", "SBR" = "Sambar deer", "PIG" = "Wild boar")

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

# Read Master Model Convergence Summary
master_df <- read.csv("Results/model_summary/Master_Model_Convergence_Summary.csv", stringsAsFactors = FALSE)

cat(sprintf("Processing covariate boxplots for %d top models...\n\n", nrow(master_df)))

active_coefs_df_list <- list()

for (i in 1:nrow(master_df)) {
  sp_name <- master_df$Species[i]
  sp_code <- master_df$Species_Code[i]
  rank <- master_df$Rank[i]
  m_type <- master_df$Type[i]
  covars_str <- master_df$Covariates[i]
  delta_val <- master_df$Delta_wAIC[i]
  
  cat(sprintf("[%d/%d] Generating covariate boxplot for %s Rank %d (%s: %s)...\n", 
              i, nrow(master_df), sp_code, rank, m_type, covars_str))
  
  rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s_Rank%d.rds", sp_code, rank)
  if (!file.exists(rds_file)) {
    rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s.rds", sp_code)
  }
  if (!file.exists(rds_file)) {
    cat("  WARNING: RDS sample file not found for", sp_code, "Rank", rank, "\n")
    next
  }
  
  s_obj <- readRDS(rds_file)
  samps <- extract_samples_matrix(s_obj)
  
  # Active covariates for this model
  active_covars <- unlist(strsplit(covars_str, " \\+ "))
  
  # Extract beta parameters
  beta_cols <- grep("^beta\\[", colnames(samps), value = TRUE)
  if (length(beta_cols) >= 7) beta_cols <- beta_cols[1:7]
  
  beta_mat <- samps[, beta_cols, drop = FALSE]
  colnames(beta_mat) <- covar_names[1:ncol(beta_mat)]
  
  # Extract inclusion indicator w parameters if present in MCMC
  w_cols <- grep("^w\\[", colnames(samps), value = TRUE)
  if (length(w_cols) >= 7) {
    w_cols <- w_cols[1:7]
    w_mat <- samps[, w_cols, drop = FALSE]
    colnames(w_mat) <- covar_names[1:ncol(w_mat)]
    w_prop <- colMeans(w_mat)
  } else {
    w_prop <- setNames(rep(0, length(covar_names)), covar_names)
    w_prop[active_covars] <- 1.0
  }
  
  # Subsample for smooth ggplot jitter plotting (600 draws)
  set.seed(42 + i)
  sub_idx <- sample(1:nrow(beta_mat), min(600, nrow(beta_mat)))
  
  # Long format data frame
  df_long <- data.frame()
  for (cov_item in covar_names) {
    if (cov_item %in% colnames(beta_mat)) {
      df_long <- rbind(df_long, data.frame(
        Covariate = cov_item,
        Value = beta_mat[sub_idx, cov_item],
        IsActive = ifelse(cov_item %in% active_covars, "Active in Model", "Unselected Prior"),
        stringsAsFactors = FALSE
      ))
    }
  }
  df_long$Covariate <- factor(df_long$Covariate, levels = covar_names)
  
  # Calculate summary stats for annotations (Mean, Median, 2.5% & 97.5% quantiles)
  summary_stats <- df_long %>%
    group_by(Covariate, IsActive) %>%
    summarise(
      Mean = mean(Value),
      Median = median(Value),
      q2.5 = quantile(Value, 0.025),
      q97.5 = quantile(Value, 0.975),
      .groups = "drop"
    )
  
  # Collect active coefficients for master comparative figure
  for (cov_item in active_covars) {
    active_coefs_df_list[[length(active_coefs_df_list) + 1]] <- data.frame(
      Species = sp_name,
      Species_Code = sp_code,
      Rank = rank,
      Covariate = cov_item,
      Value = beta_mat[sub_idx, cov_item],
      stringsAsFactors = FALSE
    )
  }
  
  # Y-axis limits expansion for clear text annotations
  y_min <- min(df_long$Value) - 0.5
  y_max <- max(df_long$Value) + 0.6
  
  # ------------------------------------------------------------------------
  # A. Image 1: Pure Covariate Beta Coefficients Boxplot with Jitter & Quantiles
  # ------------------------------------------------------------------------
  p_box <- ggplot(df_long, aes(x = Covariate, y = Value)) +
    # Boxplots
    geom_boxplot(aes(fill = IsActive), outlier.shape = NA, alpha = 0.60, width = 0.5, color = "black", linewidth = 0.6) +
    # Jittered draw points
    geom_jitter(aes(color = IsActive), width = 0.18, alpha = 0.35, size = 1.0) +
    # Reference line at y = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # 2.5% & 97.5% Quantile Error Bar Caps (Red T-bars)
    geom_errorbar(data = summary_stats, aes(x = Covariate, ymin = q2.5, ymax = q97.5),
                  inherit.aes = FALSE, width = 0.25, color = "darkred", linewidth = 0.8) +
    # Mean Mark (Red Diamond)
    geom_point(data = summary_stats, aes(x = Covariate, y = Mean),
               inherit.aes = FALSE, shape = 18, color = "red", size = 4.2) +
    # Text Annotations for 97.5% Upper Quantile
    geom_text(data = summary_stats, aes(x = Covariate, y = q97.5, label = sprintf("97.5%%: %.2f", q97.5)),
              inherit.aes = FALSE, vjust = -0.7, size = 3.2, fontface = "bold", color = "darkred") +
    # Text Annotations for 2.5% Lower Quantile
    geom_text(data = summary_stats, aes(x = Covariate, y = q2.5, label = sprintf("2.5%%: %.2f", q2.5)),
              inherit.aes = FALSE, vjust = 1.5, size = 3.2, fontface = "bold", color = "darkred") +
    # Text Annotations for Mean Value (marked with Red Diamond)
    geom_text(data = summary_stats, aes(x = Covariate, y = Mean, label = sprintf("μ = %.2f", Mean)),
              inherit.aes = FALSE, hjust = -0.25, size = 3.0, fontface = "bold", color = "red") +
    scale_fill_manual(values = c("Active in Model" = "#1F77B4", "Unselected Prior" = "#E0E0E0")) +
    scale_color_manual(values = c("Active in Model" = "#0B3C5D", "Unselected Prior" = "#777777")) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 11, color = "gray30"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
          axis.title = element_text(face = "bold"),
          legend.position = "top",
          panel.grid.major.x = element_blank()) +
    labs(title = sprintf("Covariate Beta Coefficients: %s (Rank %d)", sp_name, rank),
         subtitle = sprintf("Red Diamond = Mean (μ) | Red T-Bars & Labels = 2.5%% & 97.5%% Quantiles | Model (%s): %s (Delta wAIC = %.2f)", 
                            m_type, covars_str, delta_val),
         x = "Covariate", y = "Posterior Beta Value")
  
  png_box1 <- sprintf("Results/Covariates/Boxplot_Covariates_%s_Rank%d.png", sp_code, rank)
  png_box2 <- sprintf("Results/Covariates/Posterior_Covariates_%s_Rank%d.png", sp_code, rank)
  jpg_box1 <- sprintf("Results/Covariates/Posterior_Covariates_%s_Rank%d.jpg", sp_code, rank)
  jpg_box2 <- sprintf("Results/Covariates/Posterior_Covariates_%s_%s_Rank%d.jpg", sp_code, ifelse(m_type=="Spatial","Spatial","NoSpatial"), rank)
  
  ggsave(png_box1, p_box, width = 10, height = 6.5, dpi = 300)
  ggsave(png_box2, p_box, width = 10, height = 6.5, dpi = 300)
  ggsave(jpg_box1, p_box, width = 10, height = 6.5, dpi = 300)
  ggsave(jpg_box2, p_box, width = 10, height = 6.5, dpi = 300)
  cat(sprintf("  Saved Pure Covariate Boxplot with Jitter: %s\n", png_box1))
  
  # ------------------------------------------------------------------------
  # B. Image 2: Dedicated Inclusion Probabilities Figure
  # ------------------------------------------------------------------------
  df_w <- data.frame(
    Covariate = factor(covar_names, levels = covar_names),
    InclusionProb = as.numeric(w_prop[covar_names]),
    IsActive = ifelse(covar_names %in% active_covars, "Included (w = 1.0)", "Not Included (w = 0.0)"),
    stringsAsFactors = FALSE
  )
  
  p_incl <- ggplot(df_w, aes(x = Covariate, y = InclusionProb, fill = IsActive)) +
    geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", InclusionProb)), vjust = -0.4, fontface = "bold", size = 4.0) +
    scale_fill_manual(values = c("Included (w = 1.0)" = "#2CA02C", "Not Included (w = 0.0)" = "#E377C2")) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 11, color = "gray30"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
          axis.title = element_text(face = "bold"),
          legend.position = "top",
          panel.grid.major.x = element_blank()) +
    labs(title = sprintf("Covariate Inclusion Probabilities: %s (Rank %d)", sp_name, rank),
         subtitle = sprintf("Model (%s): %s | Active covariates in specified model formulas have fixed w = 1.0 (100%% inclusion)", m_type, covars_str),
         x = "Covariate", y = "Inclusion Probability P(w_j = 1)")
  
  png_w1 <- sprintf("Results/Covariates/Inclusion_Probabilities_%s_Rank%d.png", sp_code, rank)
  ggsave(png_w1, p_incl, width = 9, height = 5.5, dpi = 300)
  cat(sprintf("  Saved Inclusion Probabilities Plot: %s\n", png_w1))
}

# ------------------------------------------------------------------------
# C. Master Multi-Species Active Covariates Boxplot with Jitter, Mean Marks & Quantiles
# ------------------------------------------------------------------------
cat("\nGenerating Master Multi-Species Covariate Boxplot Comparison Figure...\n")

if (length(active_coefs_df_list) > 0) {
  master_coef_df <- do.call(rbind, active_coefs_df_list)
  master_rank1_df <- master_coef_df %>% filter(Rank == 1)
  
  master_summary_stats <- master_rank1_df %>%
    group_by(Species, Covariate) %>%
    summarise(
      Mean = mean(Value),
      Median = median(Value),
      q2.5 = quantile(Value, 0.025),
      q97.5 = quantile(Value, 0.975),
      .groups = "drop"
    )
  
  p_master_box <- ggplot(master_rank1_df, aes(x = Covariate, y = Value)) +
    geom_boxplot(aes(fill = Species), outlier.shape = NA, alpha = 0.60, width = 0.55, color = "black", linewidth = 0.5) +
    geom_jitter(aes(color = Species), width = 0.18, alpha = 0.25, size = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # 2.5% & 97.5% Quantile Error Bar Caps
    geom_errorbar(data = master_summary_stats, aes(x = Covariate, ymin = q2.5, ymax = q97.5),
                  inherit.aes = FALSE, width = 0.25, color = "darkred", linewidth = 0.8) +
    # Mean Mark (Red Diamond)
    geom_point(data = master_summary_stats, aes(x = Covariate, y = Mean),
               inherit.aes = FALSE, shape = 18, color = "red", size = 3.5) +
    # Text Annotations
    geom_text(data = master_summary_stats, aes(x = Covariate, y = q97.5, label = sprintf("%.2f", q97.5)),
              inherit.aes = FALSE, vjust = -0.5, size = 2.6, fontface = "bold", color = "darkred") +
    geom_text(data = master_summary_stats, aes(x = Covariate, y = q2.5, label = sprintf("%.2f", q2.5)),
              inherit.aes = FALSE, vjust = 1.3, size = 2.6, fontface = "bold", color = "darkred") +
    facet_wrap(~Species, scales = "free_x", ncol = 3) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11, color = "gray30"),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    labs(title = "Multi-Species Active Covariate Beta Coefficients (Rank 1 Models)",
         subtitle = "Red Diamond = Mean Mark (μ) | Red T-Bars & Labels = 2.5% & 97.5% Quantile Indicators",
         x = "Active Covariates", y = "Posterior Beta Value")
  
  master_png <- "Results/Covariates/Covariate_Coefficients_Boxplot_All_Species.png"
  ggsave(master_png, p_master_box, width = 14, height = 9, dpi = 300)
  cat(sprintf("  Saved Master Multi-Species Boxplot: %s\n", master_png))
}

cat("\n========================================================================\n")
cat("SUCCESS: COVARIATE-ONLY BOXPLOTS WITH MARKS & QUANTILES COMPLETED!\n")
cat("========================================================================\n")
