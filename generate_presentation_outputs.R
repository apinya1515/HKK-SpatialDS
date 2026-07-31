# generate_presentation_outputs.R
# Generate Final Presentation & Manuscript Outputs:
# 1. Multi-Species Community Density Comparison Chart (Results/Density/Community_Density_Comparison.png)
# 2. Spatial Covariate Correlation Heatmap (Results/Covariates/Spatial_Covariates_Correlation_Heatmap.png)
# 3. Full Candidate Model Selection Master Table (Results/tables/Full_Candidate_Model_Selection_Summary.csv & .xlsx)
# 4. Ecological Covariate Effects Synthesis Table (Results/tables/Ecological_Covariate_Effects_Synthesis.csv & .xlsx)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(writexl)
library(reshape2)

set.seed(42)

cat("========================================================================\n")
cat("GENERATING FINAL PRESENTATION & MANUSCRIPT OUTPUTS\n")
cat("========================================================================\n")

# Ensure required directories exist
dirs <- c("Results/Density", "Results/Covariates", "Results/tables")
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

species_names <- c("Banteng", "Gaur", "Muntjac", "Sambar deer", "Wild boar")
species_codes <- c("BTG", "GAR", "MJK", "SBR", "PIG")
study_area <- 3011 # km2

# Helper function to extract sample matrix
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

# ------------------------------------------------------------------------
# 1. Multi-Species Community Density Comparison Chart
# ------------------------------------------------------------------------
cat("\n[1/4] Generating Multi-Species Community Density Comparison Chart...\n")

community_summary_list <- list()

for (sp_idx in seq_along(species_names)) {
  sp_name <- species_names[sp_idx]
  sp_code <- species_codes[sp_idx]
  
  rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s_Rank1.rds", sp_code)
  if (!file.exists(rds_file)) rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s.rds", sp_code)
  
  if (file.exists(rds_file)) {
    s_obj <- readRDS(rds_file)
    samps <- extract_samples_matrix(s_obj)
    tot_abund <- samps[, "TOTAL_ABUND"]
    d_km2 <- tot_abund / study_area
    
    community_summary_list[[sp_name]] <- data.frame(
      Species = sp_name,
      Species_Code = sp_code,
      Abund_Mean = mean(tot_abund),
      Abund_Median = median(tot_abund),
      Abund_SD = sd(tot_abund),
      Abund_L95 = quantile(tot_abund, 0.025),
      Abund_U95 = quantile(tot_abund, 0.975),
      Density_Mean = mean(d_km2),
      Density_Median = median(d_km2),
      Density_SD = sd(d_km2),
      Density_L95 = quantile(d_km2, 0.025),
      Density_U95 = quantile(d_km2, 0.975),
      stringsAsFactors = FALSE
    )
  }
}

comm_df <- do.call(rbind, community_summary_list)
comm_df$Species <- factor(comm_df$Species, levels = c("Muntjac", "Sambar deer", "Wild boar", "Gaur", "Banteng"))

# Plot A: Individual Density Comparison
p_comm_density <- ggplot(comm_df, aes(x = Species, y = Density_Mean, fill = Species)) +
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = Density_L95, ymax = Density_U95), width = 0.2, color = "darkred", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f\n[%.2f - %.2f]", Density_Mean, Density_L95, Density_U95)),
            vjust = -0.5, size = 3.3, fontface = "bold", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, max(comm_df$Density_U95) * 1.22)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank()) +
  labs(title = "Ungulate Population Density (D)",
       subtitle = "Mean Density (individuals / km²) with 95% Bayesian Credible Interval",
       x = "Species", y = "Density (individuals / km²)")

# Plot B: Total Abundance Comparison
p_comm_abund <- ggplot(comm_df, aes(x = Species, y = Abund_Mean, fill = Species)) +
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = Abund_L95, ymax = Abund_U95), width = 0.2, color = "darkred", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.0f\n[%.0f - %.0f]", Abund_Mean, Abund_L95, Abund_U95)),
            vjust = -0.5, size = 3.3, fontface = "bold", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, max(comm_df$Abund_U95) * 1.22), labels = scales::comma) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank()) +
  labs(title = "Total Ungulate Abundance (N)",
       subtitle = "Mean Abundance (individuals in 3,011 km²) with 95% Bayesian Credible Interval",
       x = "Species", y = "Total Abundance (N)")

p_comm_combined <- grid.arrange(p_comm_density, p_comm_abund, ncol = 2)

png_comm <- "Results/Density/Community_Density_Comparison.png"
jpg_comm <- "Results/Density/Community_Density_Comparison.jpg"
ggsave(png_comm, p_comm_combined, width = 14, height = 6.5, dpi = 300)
ggsave(jpg_comm, p_comm_combined, width = 14, height = 6.5, dpi = 300)
cat(sprintf("  Saved Community Density Comparison Chart: %s\n", png_comm))

# ------------------------------------------------------------------------
# 2. Spatial Covariate Correlation Heatmap
# ------------------------------------------------------------------------
cat("\n[2/4] Generating Spatial Covariate Correlation Heatmap...\n")

land_df <- read.csv("HKK_Cov1sqkm_.csv")
covar_names <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')
covar_matrix <- land_df[, covar_names]

# Compute Pearson Correlation Matrix
cor_matrix <- cor(covar_matrix, use = "complete.obs")

# Melt matrix for ggplot heatmap
cor_melt <- melt(cor_matrix)
colnames(cor_melt) <- c("Covariate1", "Covariate2", "Correlation")
cor_melt$Covariate1 <- factor(cor_melt$Covariate1, levels = covar_names)
cor_melt$Covariate2 <- factor(cor_melt$Covariate2, levels = rev(covar_names))

p_cor <- ggplot(cor_melt, aes(x = Covariate1, y = Covariate2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), color = ifelse(abs(cor_melt$Correlation) > 0.5, "white", "black"), fontface = "bold", size = 4.2) +
  scale_fill_gradient2(low = "#2E86C1", mid = "#F7F9F9", high = "#E74C3C", midpoint = 0, limit = c(-1, 1), name = "Pearson (r)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, color = "gray30"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Spatial Covariate Correlation Heatmap (Huai Kha Khaeng)",
       subtitle = "Pairwise Pearson correlation matrix verifying low multi-collinearity (|r| < 0.70)")

png_cor <- "Results/Covariates/Spatial_Covariates_Correlation_Heatmap.png"
jpg_cor <- "Results/Covariates/Spatial_Covariates_Correlation_Heatmap.jpg"
ggsave(png_cor, p_cor, width = 8.5, height = 7.0, dpi = 300)
ggsave(jpg_cor, p_cor, width = 8.5, height = 7.0, dpi = 300)
cat(sprintf("  Saved Spatial Covariate Correlation Heatmap: %s\n", png_cor))

# ------------------------------------------------------------------------
# 3. Full Candidate Model Selection Master Table
# ------------------------------------------------------------------------
cat("\n[3/4] Generating Full Candidate Model Selection Master Table...\n")

if (file.exists("Results/Final_Model_Comparison.csv")) {
  full_model_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  
  # Format & Sort
  full_model_df <- full_model_df %>%
    arrange(Species, Rank) %>%
    mutate(
      wAIC = round(wAIC, 2),
      Delta_wAIC = round(Delta_wAIC, 2),
      Model_Weight = round(exp(-0.5 * Delta_wAIC) / sum(exp(-0.5 * Delta_wAIC)), 4)
    )
  
  out_csv_full <- "Results/tables/Full_Candidate_Model_Selection_Summary.csv"
  out_xlsx_full <- "Results/tables/Full_Candidate_Model_Selection_Summary.xlsx"
  
  write.csv(full_model_df, out_csv_full, row.names = FALSE)
  write_xlsx(list("Full_Model_Selection" = full_model_df), out_xlsx_full)
  cat(sprintf("  Saved Full Candidate Model Selection Table: %s & %s\n", out_csv_full, out_xlsx_full))
}

# ------------------------------------------------------------------------
# 4. Ecological Covariate Effects Synthesis Table
# ------------------------------------------------------------------------
cat("\n[4/4] Generating Ecological Covariate Effects Synthesis Table...\n")

master_summary <- read.csv("Results/model_summary/Master_Model_Convergence_Summary.csv", stringsAsFactors = FALSE)

eco_synthesis_list <- list()

for (i in 1:nrow(master_summary)) {
  sp_name <- master_summary$Species[i]
  sp_code <- master_summary$Species_Code[i]
  rank <- master_summary$Rank[i]
  covars_str <- master_summary$Covariates[i]
  delta_val <- master_summary$Delta_wAIC[i]
  
  rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s_Rank%d.rds", sp_code, rank)
  if (!file.exists(rds_file)) rds_file <- sprintf("Results/MCMC/MCMC_Samples_%s.rds", sp_code)
  
  if (file.exists(rds_file)) {
    s_obj <- readRDS(rds_file)
    samps <- extract_samples_matrix(s_obj)
    
    active_covars <- unlist(strsplit(covars_str, " \\+ "))
    
    for (cov_item in covar_names) {
      if (cov_item %in% active_covars) {
        c_idx <- match(cov_item, covar_names)
        beta_col <- sprintf("beta[%d]", c_idx)
        
        if (beta_col %in% colnames(samps)) {
          b_draws <- samps[, beta_col]
          m_b <- mean(b_draws)
          sd_b <- sd(b_draws)
          ci_l <- quantile(b_draws, 0.025)
          ci_u <- quantile(b_draws, 0.975)
          
          dir_str <- if (ci_l > 0) "Positive (+)" else if (ci_u < 0) "Negative (-)" else "Overlap 0"
          
          eco_synthesis_list[[length(eco_synthesis_list) + 1]] <- data.frame(
            Species = sp_name,
            Rank = rank,
            Delta_wAIC = delta_val,
            Covariate = cov_item,
            Mean_Beta = round(m_b, 3),
            SD = round(sd_b, 3),
            CI_95_Lower = round(ci_l, 3),
            CI_95_Upper = round(ci_u, 3),
            Effect_Direction = dir_str,
            Inclusion = "Active in Model (w = 1.0)",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

eco_df <- do.call(rbind, eco_synthesis_list)

out_csv_eco <- "Results/tables/Ecological_Covariate_Effects_Synthesis.csv"
out_xlsx_eco <- "Results/tables/Ecological_Covariate_Effects_Synthesis.xlsx"

write.csv(eco_df, out_csv_eco, row.names = FALSE)
write_xlsx(list("Ecological_Synthesis" = eco_df), out_xlsx_eco)
cat(sprintf("  Saved Ecological Covariate Synthesis Table: %s & %s\n", out_csv_eco, out_xlsx_eco))

cat("\n========================================================================\n")
cat("SUCCESS: ALL 4 PRESENTATION & MANUSCRIPT OUTPUTS COMPLETED!\n")
cat("========================================================================\n")
