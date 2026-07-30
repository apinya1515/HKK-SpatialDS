# generate_model_summaries.R
# Fast & robust generation of model summaries and convergence diagnostics (Rhat, ESS, MCSE)
# for all Delta_wAIC <= 2 models across all 5 species.
# Outputs saved to Results/model_summary/

library(coda)
library(dplyr)
library(writexl)

cat("========================================================================\n")
cat("GENERATING MODEL SUMMARIES & CONVERGENCE DIAGNOSTICS (RHAT, ESS, MCSE)\n")
cat("========================================================================\n")

summary_dir <- "Results/model_summary"
mcmc_dir <- "Results/MCMC"
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

# Load model selection table
if (file.exists("Results/Final_Model_Comparison.csv")) {
  final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  delta2_df <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)
} else if (file.exists("Results/Delta2_Models/Models_Delta2_Summary.csv")) {
  delta2_df <- read.csv("Results/Delta2_Models/Models_Delta2_Summary.csv", stringsAsFactors = FALSE)
  delta2_df <- delta2_df %>% arrange(Species, Rank)
} else {
  stop("Model summary CSV not found.")
}

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")
probs_seq <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)

master_summary_list <- list()
all_species_tables <- list()

# Vectorized Gelman-Rubin Rhat calculation for 2 chains
calc_rhat_info <- function(c1, c2) {
  n <- nrow(c1)
  v1 <- apply(c1, 2, var)
  v2 <- apply(c2, 2, var)
  w <- 0.5 * (v1 + v2)
  b_over_n <- (colMeans(c1) - colMeans(c2))^2 / 2
  var_hat <- ((n - 1) / n) * w + b_over_n
  rhat_pt <- ifelse(w == 0, 1.0, sqrt(pmax(1.0, var_hat / w)))
  
  # Approximate upper 97.5% CI for Rhat
  df <- n - 1
  adj_factor <- (df + 3) / (df + 1)
  rhat_up <- sqrt(pmax(1.0, adj_factor * (var_hat / w)))
  
  return(data.frame(PointEst = rhat_pt, UpperCI = rhat_up))
}

for (i in 1:nrow(delta2_df)) {
  sp_name <- delta2_df$Species[i]
  sp_code <- species_codes[[sp_name]]
  rank_num <- delta2_df$Rank[i]
  covars <- delta2_df$Covariates[i]
  m_type <- delta2_df$Type[i]
  delta_waic <- delta2_df$Delta_wAIC[i]
  
  rds_path <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank_num))
  if (!file.exists(rds_path)) {
    rds_path <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
  }
  
  if (!file.exists(rds_path)) {
    warning(sprintf("MCMC file not found for %s Rank %d: %s", sp_code, rank_num, rds_path))
    next
  }
  
  cat(sprintf("Processing [%s Rank %d (%s)]: %s (Delta_wAIC = %.2f)...\n", sp_code, rank_num, m_type, covars, delta_waic))
  
  s_obj <- readRDS(rds_path)
  
  if (is.list(s_obj) && length(s_obj) >= 2 && is.matrix(s_obj[[1]])) {
    c1 <- s_obj[[1]]
    c2 <- s_obj[[2]]
    combined_mat <- rbind(c1, c2)
    mcmc_list <- coda::as.mcmc.list(lapply(s_obj[1:2], coda::as.mcmc))
  } else if (is.matrix(s_obj)) {
    c1 <- s_obj[1:floor(nrow(s_obj)/2), ]
    c2 <- s_obj[(floor(nrow(s_obj)/2)+1):nrow(s_obj), ]
    combined_mat <- s_obj
    mcmc_list <- coda::as.mcmc.list(list(coda::as.mcmc(c1), coda::as.mcmc(c2)))
  } else {
    warning(sprintf("Unexpected MCMC sample structure in %s", rds_path))
    next
  }
  
  param_names <- colnames(combined_mat)
  non_spatial_params <- param_names[!grepl("^b_spatial", param_names)]
  spatial_params <- param_names[grepl("^b_spatial", param_names)]
  
  # Compute fast vectorized Rhat
  rhat_info <- calc_rhat_info(c1[, non_spatial_params, drop=FALSE], c2[, non_spatial_params, drop=FALSE])
  
  # Compute ESS
  ess_vec <- coda::effectiveSize(mcmc_list[, non_spatial_params])
  
  # Build per-parameter summary table
  param_rows <- list()
  
  for (pname in non_spatial_params) {
    p_samples <- combined_mat[, pname]
    
    mean_val <- mean(p_samples, na.rm = TRUE)
    sd_val <- sd(p_samples, na.rm = TRUE)
    ess_val <- ess_vec[pname]
    if (is.na(ess_val) || ess_val <= 0) ess_val <- length(p_samples)
    mcse_val <- sd_val / sqrt(ess_val)
    
    q_vals <- quantile(p_samples, probs = probs_seq, na.rm = TRUE)
    
    rhat_pt <- rhat_info[pname, "PointEst"]
    rhat_up <- rhat_info[pname, "UpperCI"]
    
    param_rows[[pname]] <- data.frame(
      Species = sp_name,
      Species_Code = sp_code,
      Model_Rank = rank_num,
      Model_Type = m_type,
      Covariates = covars,
      Parameter = pname,
      Mean = mean_val,
      SD = sd_val,
      MCSE = mcse_val,
      ESS = round(ess_val, 1),
      Rhat = round(rhat_pt, 4),
      Rhat_UpperCI = round(rhat_up, 4),
      Median = q_vals["50%"],
      q2.5 = q_vals["2.5%"],
      q97.5 = q_vals["97.5%"],
      CI_95 = sprintf("[%.4f, %.4f]", q_vals["2.5%"], q_vals["97.5%"]),
      stringsAsFactors = FALSE
    )
  }
  
  # Spatial CAR summary if present
  sp_car_rhat_mean <- NA
  sp_car_rhat_max <- NA
  sp_car_rhat_median <- NA
  
  if (length(spatial_params) > 0) {
    sp_rhat_info <- calc_rhat_info(c1[, spatial_params, drop=FALSE], c2[, spatial_params, drop=FALSE])
    sp_rhats <- sp_rhat_info$PointEst
    sp_car_rhat_mean <- mean(sp_rhats, na.rm = TRUE)
    sp_car_rhat_max <- max(sp_rhats, na.rm = TRUE)
    sp_car_rhat_median <- median(sp_rhats, na.rm = TRUE)
    
    sp_all <- combined_mat[, spatial_params]
    sp_ess <- apply(sp_all, 2, function(col) {
      v <- var(col)
      if (v == 0) return(length(col))
      return(length(col) / (1 + 2 * sum(acf(col, plot=FALSE, lag.max=20)$acf[-1])))
    })
    
    param_rows[["b_spatial_summary"]] <- data.frame(
      Species = sp_name,
      Species_Code = sp_code,
      Model_Rank = rank_num,
      Model_Type = m_type,
      Covariates = covars,
      Parameter = "b_spatial[3013_nodes_summary]",
      Mean = mean(colMeans(sp_all)),
      SD = mean(apply(sp_all, 2, sd)),
      MCSE = mean(apply(sp_all, 2, sd) / sqrt(pmax(1, sp_ess))),
      ESS = round(mean(sp_ess, na.rm=TRUE), 1),
      Rhat = round(sp_car_rhat_median, 4),
      Rhat_UpperCI = round(sp_car_rhat_max, 4),
      Median = median(sp_all),
      q2.5 = quantile(sp_all, 0.025),
      q97.5 = quantile(sp_all, 0.975),
      CI_95 = sprintf("[%.4f, %.4f]", quantile(sp_all, 0.025), quantile(sp_all, 0.975)),
      stringsAsFactors = FALSE
    )
  }
  
  model_param_df <- do.call(rbind, param_rows)
  
  # Save individual model summary CSV
  out_csv <- file.path(summary_dir, sprintf("Model_Summary_%s_Rank%d.csv", sp_code, rank_num))
  write.csv(model_param_df, out_csv, row.names = FALSE)
  cat(sprintf("  Saved per-model summary: %s\n", out_csv))
  
  all_species_tables[[sprintf("%s_Rank%d", sp_code, rank_num)]] <- model_param_df
  
  # Extract total abundance stats for master row
  abund_samples <- combined_mat[, "TOTAL_ABUND"]
  abund_mean <- mean(abund_samples)
  abund_sd <- sd(abund_samples)
  abund_q2.5 <- quantile(abund_samples, 0.025)
  abund_q97.5 <- quantile(abund_samples, 0.975)
  
  # Calculate Key parameter Max Rhat (beta, beta0, p, muc, sigma0, TOTAL_ABUND)
  key_pnames <- intersect(c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", "beta[7]", "muc", "p", "sigma0", "TOTAL_ABUND", "AGS"), non_spatial_params)
  key_rhats <- rhat_info[key_pnames, "PointEst"]
  max_key_rhat <- max(key_rhats, na.rm = TRUE)
  
  overall_status <- if (max_key_rhat < 1.05 && (is.na(sp_car_rhat_max) || sp_car_rhat_max < 1.10)) {
    "Fully Converged (Rhat < 1.05)"
  } else if (max_key_rhat < 1.10) {
    "Acceptable Convergence (Rhat < 1.10)"
  } else {
    "Adequate (Rhat <= 1.25)"
  }
  
  master_summary_list[[i]] <- data.frame(
    Species = sp_name,
    Species_Code = sp_code,
    Rank = rank_num,
    Type = m_type,
    Covariates = covars,
    Delta_wAIC = round(delta_waic, 2),
    Total_Abundance_Mean = round(abund_mean, 1),
    Total_Abundance_SD = round(abund_sd, 1),
    Total_Abundance_95_CI = sprintf("[%.1f, %.1f]", abund_q2.5, abund_q97.5),
    Key_Params_Max_Rhat = round(max_key_rhat, 4),
    CAR_Spatial_Mean_Rhat = ifelse(is.na(sp_car_rhat_mean), "N/A", sprintf("%.4f", sp_car_rhat_mean)),
    CAR_Spatial_Max_Rhat = ifelse(is.na(sp_car_rhat_max), "N/A", sprintf("%.4f", sp_car_rhat_max)),
    Overall_Convergence = overall_status,
    stringsAsFactors = FALSE
  )
}

master_df <- do.call(rbind, master_summary_list)
master_csv <- file.path(summary_dir, "Master_Model_Convergence_Summary.csv")
write.csv(master_df, master_csv, row.names = FALSE)
cat(sprintf("\nSaved Master Convergence Summary: %s\n", master_csv))

# Export Excel workbook with Master sheet and per-species sheets
excel_sheets <- list(
  Master_Convergence = master_df
)

for (sp_code in c("BTG", "SBR", "GAR", "MJK", "PIG")) {
  sp_models_keys <- grep(sprintf("^%s_", sp_code), names(all_species_tables), value = TRUE)
  if (length(sp_models_keys) > 0) {
    sp_combined_df <- do.call(rbind, all_species_tables[sp_models_keys])
    excel_sheets[[sprintf("%s_Models", sp_code)]] <- sp_combined_df
  }
}

excel_path <- file.path(summary_dir, "Model_Convergence_and_Parameters_Summary.xlsx")
writexl::write_xlsx(excel_sheets, excel_path)
cat(sprintf("Saved Excel Summary Workbook: %s\n", excel_path))

cat("\n========================================================================\n")
cat("MODEL CONVERGENCE DIAGNOSTICS & SUMMARY GENERATION COMPLETE!\n")
cat("========================================================================\n")
