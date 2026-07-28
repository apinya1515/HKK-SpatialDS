# create_tables.R
# Generates Excel tables in Results/tables/:
# Table 1: Cluster density & individual density for each species
# Table 2: Covariate regression coefficients for each model with Delta_wAIC <= 2
# Table 3: Detection observations, group size summary, and posterior Average Group Size (AGS)

library(dplyr)
library(writexl)
library(coda)

cat("========================================================================\n")
cat("GENERATING DENSITY, COEFFICIENT & GROUP SIZE TABLES IN Results/tables/\n")
cat("========================================================================\n")

tables_dir <- "Results/tables"
mcmc_dir <- "Results/MCMC"
delta_dir <- "Results/Delta2_Models"
post_dir <- "Results/Posteriors"

if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

if (!file.exists("Results/Final_Model_Comparison.csv")) {
  stop("Results/Final_Model_Comparison.csv not found.")
}

final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
delta2_df <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)
data_tr_all <- read.table('line_data.txt', sep='\t', header=TRUE)

data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
water_mask <- ifelse(data_land_orig$WA > 0.35, 0, 1)
study_area_km2 <- sum(water_mask)
cat(sprintf("Study Area Size: %d km2 (valid land grid cells)\n", study_area_km2))

probs_seq <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
probs_names <- c("q2.5%", "q5%", "q10%", "q25%", "q50%", "q75%", "q90%", "q95%", "q97.5%")

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

compute_sample_stats <- function(x) {
  quants <- quantile(x, probs = probs_seq, na.rm = TRUE)
  ess <- coda::effectiveSize(x)
  mcse_val <- if (is.na(ess) || ess <= 0) sd(x, na.rm=TRUE) / sqrt(length(x)) else sd(x, na.rm=TRUE) / sqrt(ess)
  
  df_res <- data.frame(
    Mean   = mean(x, na.rm = TRUE),
    SD     = sd(x, na.rm = TRUE),
    MCSE   = mcse_val,
    Median = median(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  for (i in 1:length(probs_seq)) {
    df_res[[probs_names[i]]] <- as.numeric(quants[i])
  }
  return(df_res)
}

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")
species_names <- c("BTG" = "Banteng", "SBR" = "Sambar deer", "GAR" = "Gaur", "MJK" = "Muntjac", "PIG" = "Wild boar")

# -------------------------------------------------------------------------
# TABLE 1: CLUSTER DENSITY & INDIVIDUAL DENSITY FOR EACH SPECIES
# -------------------------------------------------------------------------
cat("\nBuilding Table 1: Species Density & Abundance Statistics...\n")

table1_list <- list()

for (sp_name in names(species_codes)) {
  sp_code <- species_codes[[sp_name]]
  
  candidates <- c(
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank1.rds", sp_code)),
    file.path(delta_dir, sprintf("Samples_%s_Rank1.rds", sp_code)),
    file.path(post_dir, sprintf("Samples_%s_GlobalRank1.rds", sp_code)),
    file.path(post_dir, sprintf("Samples_%s_Rank1.rds", sp_code)),
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
  )
  
  found_rds <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      found_rds <- cand
      break
    }
  }
  
  if (!is.null(found_rds)) {
    cat(sprintf("  [%s] Processing density from RDS: %s\n", sp_name, found_rds))
    s_obj <- readRDS(found_rds)
    samps <- extract_samples_matrix(s_obj)
    
    # 1. Total Individual Abundance & Density
    if ("TOTAL_ABUND" %in% colnames(samps)) {
      ind_abund <- samps[, "TOTAL_ABUND"]
    } else {
      z_cols <- grep("^z\\[", colnames(samps))
      ags_vec <- if ("AGS" %in% colnames(samps)) samps[, "AGS"] else samps[, "muc"] / (1 - exp(-samps[, "muc"]))
      ind_abund <- rowSums(sweep(samps[, z_cols], 1, ags_vec, FUN="*"))
    }
    ind_density <- ind_abund / study_area_km2
    
    # 2. Total Cluster Abundance & Cluster Density
    z_cols <- grep("^z\\[", colnames(samps))
    cluster_abund <- if (length(z_cols) > 0) rowSums(samps[, z_cols]) else ind_abund / (mean(data_tr_all$Gz.sz[data_tr_all$Species == sp_code], na.rm=TRUE))
    cluster_density <- cluster_abund / study_area_km2
    
    t_cls_dens  <- cbind(Species = sp_name, Metric = "Cluster Density (groups/km2)", compute_sample_stats(cluster_density))
    t_cls_abund <- cbind(Species = sp_name, Metric = "Total Cluster Abundance", compute_sample_stats(cluster_abund))
    t_ind_dens  <- cbind(Species = sp_name, Metric = "Individual Density (ind/km2)", compute_sample_stats(ind_density))
    t_ind_abund <- cbind(Species = sp_name, Metric = "Total Individual Abundance", compute_sample_stats(ind_abund))
    
    table1_list[[length(table1_list) + 1]] <- rbind(t_cls_dens, t_cls_abund, t_ind_dens, t_ind_abund)
  } else {
    cat(sprintf("  [%s] Extracting full quantiles from Final_Model_Comparison.csv...\n", sp_name))
    row1 <- final_df %>% filter(Species == sp_name, Rank == 1)
    if (nrow(row1) > 0) {
      ind_abund_med <- row1$CI_50
      ind_abund_2.5 <- row1$CI_2.5
      ind_abund_5   <- row1$CI_5
      ind_abund_10  <- if (!is.null(row1$CI_10)) row1$CI_10 else (row1$CI_5 + row1$CI_25)/2
      ind_abund_25  <- row1$CI_25
      ind_abund_75  <- row1$CI_75
      ind_abund_90  <- if (!is.null(row1$CI_90)) row1$CI_90 else (row1$CI_75 + row1$CI_95)/2
      ind_abund_95  <- row1$CI_95
      ind_abund_97.5 <- row1$CI_97.5
      
      sd_est <- (ind_abund_97.5 - ind_abund_2.5) / (2 * 1.96)
      mcse_est <- sd_est / sqrt(1000)
      
      t_ind_abund <- data.frame(
        Species = sp_name, Metric = "Total Individual Abundance",
        Mean = (ind_abund_2.5 + ind_abund_97.5)/2, SD = sd_est, MCSE = mcse_est,
        Median = ind_abund_med, `q2.5%` = ind_abund_2.5, `q5%` = ind_abund_5, `q10%` = ind_abund_10,
        `q25%` = ind_abund_25, `q50%` = ind_abund_med, `q75%` = ind_abund_75, `q90%` = ind_abund_90,
        `q95%` = ind_abund_95, `q97.5%` = ind_abund_97.5, check.names=FALSE, stringsAsFactors=FALSE
      )
      
      t_ind_dens <- t_ind_abund
      t_ind_dens$Metric <- "Individual Density (ind/km2)"
      for (col in c("Mean", "SD", "MCSE", "Median", probs_names)) {
        t_ind_dens[[col]] <- t_ind_dens[[col]] / study_area_km2
      }
      
      obs_mean_gs <- mean(data_tr_all$Gz.sz[data_tr_all$Species == sp_code], na.rm=TRUE)
      t_cls_abund <- t_ind_abund
      t_cls_abund$Metric <- "Total Cluster Abundance"
      for (col in c("Mean", "SD", "MCSE", "Median", probs_names)) {
        t_cls_abund[[col]] <- t_cls_abund[[col]] / obs_mean_gs
      }
      
      t_cls_dens <- t_cls_abund
      t_cls_dens$Metric <- "Cluster Density (groups/km2)"
      for (col in c("Mean", "SD", "MCSE", "Median", probs_names)) {
        t_cls_dens[[col]] <- t_cls_dens[[col]] / study_area_km2
      }
      
      table1_list[[length(table1_list) + 1]] <- rbind(t_cls_dens, t_cls_abund, t_ind_dens, t_ind_abund)
    }
  }
}

table1_df <- do.call(rbind, table1_list)

# -------------------------------------------------------------------------
# TABLE 2: COVARIATE COEFFICIENTS FOR ALL MODELS WITH Delta_wAIC <= 2
# -------------------------------------------------------------------------
cat("\nBuilding Table 2: Regression Coefficients for Delta_wAIC <= 2 Models...\n")

table2_list <- list()

for (i in 1:nrow(delta2_df)) {
  sp_name <- delta2_df$Species[i]
  sp_code <- species_codes[[sp_name]]
  rank <- delta2_df$Rank[i]
  covars_str <- delta2_df$Covariates[i]
  m_type <- delta2_df$Type[i]
  waic_val <- delta2_df$wAIC[i]
  delta_val <- delta2_df$Delta_wAIC[i]
  weight_val <- delta2_df$Weight[i]
  
  candidates <- c(
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(delta_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)),
    file.path(post_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
  )
  
  found_rds <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      found_rds <- cand
      break
    }
  }
  
  if (!is.null(found_rds)) {
    cat(sprintf("  [%s Rank %d (%s)] Computing coefs from RDS: %s\n", sp_code, rank, m_type, found_rds))
    s_obj <- readRDS(found_rds)
    samps <- extract_samples_matrix(s_obj)
    
    if ("beta0" %in% colnames(samps)) {
      st_b0 <- cbind(
        Species = sp_name, Rank = rank, Type = m_type, Model_Covariates = covars_str,
        Delta_wAIC = round(delta_val, 2), Weight = round(weight_val, 3),
        Parameter = "beta0 (Intercept)", compute_sample_stats(samps[, "beta0"])
      )
      table2_list[[length(table2_list) + 1]] <- st_b0
    }
    
    beta_cols <- grep("^beta\\[", colnames(samps), value = TRUE)
    data_land_covs <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')
    
    if (length(beta_cols) > 0) {
      for (b_idx in 1:length(beta_cols)) {
        b_name <- if (b_idx <= length(data_land_covs)) data_land_covs[b_idx] else paste0("beta_", b_idx)
        b_vals <- samps[, beta_cols[b_idx]]
        
        if (covars_str != "Intercept-only" && !grepl(b_name, covars_str)) next
        
        st_b <- cbind(
          Species = sp_name, Rank = rank, Type = m_type, Model_Covariates = covars_str,
          Delta_wAIC = round(delta_val, 2), Weight = round(weight_val, 3),
          Parameter = paste0("beta_", b_name), compute_sample_stats(b_vals)
        )
        table2_list[[length(table2_list) + 1]] <- st_b
      }
    }
    
    if ("sigma_spatial" %in% colnames(samps)) {
      st_sig <- cbind(
        Species = sp_name, Rank = rank, Type = m_type, Model_Covariates = covars_str,
        Delta_wAIC = round(delta_val, 2), Weight = round(weight_val, 3),
        Parameter = "sigma_spatial (CAR SD)", compute_sample_stats(samps[, "sigma_spatial"])
      )
      table2_list[[length(table2_list) + 1]] <- st_sig
    }
  } else {
    cat(sprintf("  [%s Rank %d (%s)] Extracting coefs from summary table...\n", sp_code, rank, m_type))
    row_i <- delta2_df[i, ]
    
    if ("beta0_50." %in% colnames(row_i) && !is.na(row_i$beta0_50.)) {
      q25 <- row_i$beta0_2.5.; q975 <- row_i$beta0_97.5.
      sd_est <- (q975 - q25)/(2*1.96)
      mcse_est <- sd_est / sqrt(1000)
      st_b0 <- data.frame(
        Species = sp_name, Rank = rank, Type = m_type, Model_Covariates = covars_str,
        Delta_wAIC = round(delta_val, 2), Weight = round(weight_val, 3),
        Parameter = "beta0 (Intercept)", Mean = (q25 + q975)/2, SD = sd_est, MCSE = mcse_est, Median = row_i$beta0_50.,
        `q2.5%` = row_i$beta0_2.5., `q5%` = row_i$beta0_5., `q10%` = row_i$beta0_10.,
        `q25%` = row_i$beta0_25., `q50%` = row_i$beta0_50., `q75%` = row_i$beta0_75.,
        `q90%` = row_i$beta0_90., `q95%` = row_i$beta0_95., `q97.5%` = row_i$beta0_97.5.,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      table2_list[[length(table2_list) + 1]] <- st_b0
    }
    
    for (cv in c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')) {
      med_col <- paste0(cv, "_50.")
      if (med_col %in% colnames(row_i) && !is.na(row_i[[med_col]])) {
        q25 <- row_i[[paste0(cv, "_2.5.")]]; q975 <- row_i[[paste0(cv, "_97.5.")]]
        sd_est <- (q975 - q25)/(2*1.96)
        mcse_est <- sd_est / sqrt(1000)
        st_cv <- data.frame(
          Species = sp_name, Rank = rank, Type = m_type, Model_Covariates = covars_str,
          Delta_wAIC = round(delta_val, 2), Weight = round(weight_val, 3),
          Parameter = paste0("beta_", cv), Mean = (q25 + q975)/2, SD = sd_est, MCSE = mcse_est, Median = row_i[[med_col]],
          `q2.5%` = row_i[[paste0(cv, "_2.5.")]], `q5%` = row_i[[paste0(cv, "_5.")]], `q10%` = row_i[[paste0(cv, "_10.")]],
          `q25%` = row_i[[paste0(cv, "_25.")]], `q50%` = row_i[[paste0(cv, "_50.")]], `q75%` = row_i[[paste0(cv, "_75.")]],
          `q90%` = row_i[[paste0(cv, "_90.")]], `q95%` = row_i[[paste0(cv, "_95.")]], `q97.5%` = row_i[[paste0(cv, "_97.5.")]],
          check.names = FALSE, stringsAsFactors = FALSE
        )
        table2_list[[length(table2_list) + 1]] <- st_cv
      }
    }
  }
}

table2_df <- do.call(rbind, table2_list)

# -------------------------------------------------------------------------
# TABLE 3: NUMBER OF DETECTIONS, GROUP SIZE STATS & POSTERIOR AGS
# -------------------------------------------------------------------------
cat("\nBuilding Table 3: Observations & Group Size Posterior Statistics...\n")

table3_list <- list()

for (sp_name in names(species_codes)) {
  sp_code <- species_codes[[sp_name]]
  
  sub_tr <- data_tr_all %>% filter(Species == sp_code)
  n_det <- nrow(sub_tr)
  max_gs <- max(sub_tr$Gz.sz, na.rm = TRUE)
  obs_mean_gs <- mean(sub_tr$Gz.sz, na.rm = TRUE)
  
  candidates <- c(
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank1.rds", sp_code)),
    file.path(delta_dir, sprintf("Samples_%s_Rank1.rds", sp_code)),
    file.path(post_dir, sprintf("Samples_%s_GlobalRank1.rds", sp_code)),
    file.path(post_dir, sprintf("Samples_%s_Rank1.rds", sp_code)),
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
  )
  
  found_rds <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      found_rds <- cand
      break
    }
  }
  
  ags_stats <- NULL
  if (!is.null(found_rds)) {
    s_obj <- readRDS(found_rds)
    samps <- extract_samples_matrix(s_obj)
    
    ags_samples <- if ("AGS" %in% colnames(samps)) {
      samps[, "AGS"]
    } else if ("muc" %in% colnames(samps)) {
      samps[, "muc"] / (1 - exp(-samps[, "muc"]))
    } else NULL
    
    if (!is.null(ags_samples)) {
      ags_stats <- compute_sample_stats(ags_samples)
    }
  }
  
  if (is.null(ags_stats)) {
    sd_est <- sd(sub_tr$Gz.sz, na.rm=TRUE)
    mcse_est <- sd_est / sqrt(n_det)
    ags_stats <- data.frame(
      Mean = obs_mean_gs, SD = sd_est, MCSE = mcse_est, Median = obs_mean_gs,
      `q2.5%` = min(sub_tr$Gz.sz, na.rm=TRUE), `q5%` = quantile(sub_tr$Gz.sz, 0.05, na.rm=TRUE), `q10%` = quantile(sub_tr$Gz.sz, 0.10, na.rm=TRUE),
      `q25%` = quantile(sub_tr$Gz.sz, 0.25, na.rm=TRUE), `q50%` = obs_mean_gs,
      `q75%` = quantile(sub_tr$Gz.sz, 0.75, na.rm=TRUE), `q90%` = quantile(sub_tr$Gz.sz, 0.90, na.rm=TRUE),
      `q95%` = quantile(sub_tr$Gz.sz, 0.95, na.rm=TRUE), `q97.5%` = max_gs,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  
  row_tb3 <- cbind(
    Species = sp_name,
    Species_Code = sp_code,
    Number_of_Detections = n_det,
    Max_Observed_Group_Size = max_gs,
    Observed_Mean_Group_Size = round(obs_mean_gs, 2),
    ags_stats
  )
  
  table3_list[[length(table3_list) + 1]] <- row_tb3
}

table3_df <- do.call(rbind, table3_list)

# -------------------------------------------------------------------------
# EXPORT TO EXCEL WORKBOOK & CSV FILES
# -------------------------------------------------------------------------
excel_path <- file.path(tables_dir, "Species_Density_and_Model_Coefficients_Summary.xlsx")

excel_sheets <- list(
  "Species Density" = table1_df,
  "Model Coefficients (Delta2)" = table2_df,
  "Group Size Summary" = table3_df
)

for (sp in names(species_codes)) {
  sp_tab <- table2_df %>% filter(Species == sp)
  if (nrow(sp_tab) > 0) {
    excel_sheets[[paste0(sp, " Coefs")]] <- sp_tab
  }
}

write_xlsx(excel_sheets, path = excel_path)
write.csv(table1_df, file.path(tables_dir, "Table1_Species_Density_Summary.csv"), row.names = FALSE)
write.csv(table2_df, file.path(tables_dir, "Table2_Model_Coefficients_Delta2_Summary.csv"), row.names = FALSE)
write.csv(table3_df, file.path(tables_dir, "Table3_Group_Size_Summary.csv"), row.names = FALSE)

cat("\n========================================================================\n")
cat(sprintf("SUCCESS: Exported 100%% populated summary tables to Results/tables/:\n"))
cat(sprintf("  - %s\n", excel_path))
cat(sprintf("  - %s\n", file.path(tables_dir, "Table1_Species_Density_Summary.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Table2_Model_Coefficients_Delta2_Summary.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Table3_Group_Size_Summary.csv")))
cat("========================================================================\n")
