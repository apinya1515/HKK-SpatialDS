# generate_detection_outputs.R
# Creates Results/detection/ and generates:
# 1. Diagram plots of detection for EVERY model with Delta_wAIC <= 2 across ALL 5 species
# 2. Excel file (and CSV) of all related parameters of detection function with:
#    Mean, SD, Monte Carlo SE, Median, and 95% HPD Interval
# 3. Multi-species comparative figure: Detection_Plot_All_Species_Comparison.png
# 4. Cleans up any duplicate filename aliases (.jpg or unranked duplicates)

library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)
library(writexl)

set.seed(42)

detection_dir <- "Results/detection"
mcmc_dir <- "Results/MCMC"
post_dir <- "Results/Posteriors"

if (!dir.exists(detection_dir)) dir.create(detection_dir, recursive = TRUE)
if (!dir.exists(mcmc_dir)) dir.create(mcmc_dir, recursive = TRUE)

data_tr_all <- read.table('line_data.txt', sep='\t', header=TRUE)

species_info <- list(
  BTG = "Banteng",
  SBR = "Sambar deer",
  GAR = "Gaur",
  MJK = "Muntjac",
  PIG = "Wild boar"
)

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

dist_limit <- 100

calc_mcse <- function(x) {
  ess <- coda::effectiveSize(x)
  if (is.na(ess) || ess <= 0) return(sd(x, na.rm=TRUE) / sqrt(length(x)))
  return(sd(x, na.rm=TRUE) / sqrt(ess))
}

calc_hpd <- function(x) {
  tryCatch({
    hpd <- coda::HPDinterval(coda::as.mcmc(x), prob = 0.95)
    c(Lower = as.numeric(hpd[1, "lower"]), Upper = as.numeric(hpd[1, "upper"]))
  }, error = function(e) {
    c(Lower = as.numeric(quantile(x, 0.025, na.rm=TRUE)), Upper = as.numeric(quantile(x, 0.975, na.rm=TRUE)))
  })
}

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

# Read model selection results from Final_Model_Comparison.csv
if (file.exists("Results/Final_Model_Comparison.csv")) {
  final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  delta2_models <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)
} else {
  stop("Results/Final_Model_Comparison.csv not found.")
}

cat("========================================================================\n")
cat(sprintf("GENERATING CLEAN DETECTION PLOTS FOR ALL %d DELTA_wAIC <= 2 MODELS\n", nrow(delta2_models)))
cat("========================================================================\n")

# Clean existing directory to ensure no duplicate alias files remain
existing_files <- list.files(detection_dir, full.names = TRUE)
if (length(existing_files) > 0) unlink(existing_files)

all_summary_list <- list()
rank1_summary_list <- list()

for (i in 1:nrow(delta2_models)) {
  sp_name <- delta2_models$Species[i]
  sp_code <- species_codes[[sp_name]]
  rank <- delta2_models$Rank[i]
  m_type <- delta2_models$Type[i]
  covars_str <- delta2_models$Covariates[i]
  delta_val <- delta2_models$Delta_wAIC[i]
  
  cat(sprintf("\n--- Processing Detection for [%s Rank %d (%s)]: %s (Delta_wAIC = %.2f) ---\n", 
              sp_code, rank, m_type, covars_str, delta_val))
  
  candidates <- c(
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(post_dir, sprintf("Samples_%s_GlobalRank%d.rds", sp_code, rank)),
    file.path(post_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank)),
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code)),
    file.path(mcmc_dir, sprintf("MCMC_Samples_%s_NoSpatial.rds", sp_code))
  )
  
  found_rds <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      found_rds <- cand
      break
    }
  }
  
  samples_matrix <- NULL
  if (!is.null(found_rds)) {
    cat(sprintf("  Loading MCMC samples from: %s\n", found_rds))
    samples_matrix <- tryCatch({
      s_obj <- readRDS(found_rds)
      extract_samples_matrix(s_obj)
    }, error = function(e) {
      cat(sprintf("  Warning reading %s: %s. Falling back to Rank 1 RDS...\n", found_rds, e$message))
      fallback_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank1.rds", sp_code))
      if (!file.exists(fallback_rds)) fallback_rds <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s.rds", sp_code))
      if (file.exists(fallback_rds)) {
        extract_samples_matrix(readRDS(fallback_rds))
      } else NULL
    })
  }
  
  if (is.null(samples_matrix)) {
    cat(sprintf("  WARNING: MCMC sample matrix for %s Rank %d not found. Skipping plot...\n", sp_code, rank))
    next
  }
  
  det_param_cols <- grep("^(sigma0|p|sigma\\[|pi\\[|gs_k\\[|muc)", colnames(samples_matrix), value = TRUE)
  
  if (length(det_param_cols) > 0) {
    sp_summaries <- list()
    for (param in det_param_cols) {
      vals <- samples_matrix[, param]
      hpd <- calc_hpd(vals)
      mcse_val <- calc_mcse(vals)
      
      gs_cat <- "-"
      dist_cat <- "-"
      if (grepl("sigma\\[", param) || grepl("gs_k\\[", param)) {
        k_idx <- gsub(".*\\[([0-9]+)\\].*", "\\1", param)
        gs_cat <- paste("Group Size Class", k_idx)
      } else if (grepl("pi\\[", param)) {
        matches <- regmatches(param, regexec("pi\\[([0-9]+),\\s*([0-9]+)\\]", param))[[1]]
        if (length(matches) == 3) {
          gs_cat <- paste("Group Size Class", matches[2])
          dist_cat <- paste("Distance Bin", matches[3])
        }
      }
      
      sp_summaries[[length(sp_summaries) + 1]] <- data.frame(
        Species = sp_name,
        Species_Code = sp_code,
        Rank = rank,
        Type = m_type,
        Covariates = covars_str,
        Parameter = param,
        Group_Size_Class = gs_cat,
        Distance_Bin = dist_cat,
        Mean = mean(vals),
        SD = sd(vals),
        MCSE = mcse_val,
        Median = median(vals),
        HPD_95_Lower = hpd["Lower"],
        HPD_95_Upper = hpd["Upper"],
        stringsAsFactors = FALSE
      )
    }
    
    sp_summary_df <- do.call(rbind, sp_summaries)
    all_summary_list[[length(all_summary_list) + 1]] <- sp_summary_df
    if (rank == 1) rank1_summary_list[[sp_code]] <- sp_summary_df
  }
  
  sigma_cols <- grep("^sigma\\[", colnames(samples_matrix), value = TRUE)
  if (length(sigma_cols) > 0) {
    sigma_indices <- as.numeric(gsub("sigma\\[([0-9]+)\\]", "\\1", sigma_cols))
    sigma_cols <- sigma_cols[order(sigma_indices)]
    
    K_classes <- length(sigma_cols)
    sigma_medians <- apply(samples_matrix[, sigma_cols, drop = FALSE], 2, median)
    # Check monotonicity: sigma should increase with group size class (larger groups detected farther)
    if (any(diff(sigma_medians) < 0)) {
      cat(sprintf("  WARNING: sigma values are NOT monotonically increasing for %s Rank %d!\n", sp_code, rank))
      cat(sprintf("    sigma medians: %s\n", paste(round(sigma_medians, 2), collapse = ", ")))
      cat("    This may indicate the p parameter posterior includes negative values.\n")
    }
    
    gs_labels <- paste("Group Size Class", 1:K_classes)
    colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")
    if (K_classes != 5) colors <- rainbow(K_classes, s = 0.8, v = 0.7)
    
    # Save clean unique plot file: Detection_Plot_[Species_Code]_Rank[N].png
    plot_png_rank <- file.path(detection_dir, sprintf("Detection_Plot_%s_Rank%d.png", sp_code, rank))
    
    png(plot_png_rank, width = 2400, height = 1800, res = 300)
    x_seq <- seq(0, dist_limit, length.out = 200)
    
    plot(NULL, xlim = c(0, dist_limit), ylim = c(0, 1.05),
         xlab = "Distance from Transect Line (m)", ylab = "Detection Probability g(x)",
         main = sprintf("Detection Function: %s (Rank %d - %s)", sp_name, rank, m_type),
         sub = sprintf("Covariates: %s | Delta_wAIC = %.2f", covars_str, delta_val),
         cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3, cex.sub = 1.0, font.main = 2,
         las = 1, bty = "l")
    
    grid(col = "gray85", lty = "dotted")
    
    for (k in 1:K_classes) {
      sig <- sigma_medians[k]
      gx <- exp(- (x_seq^2) / (2 * sig^2))
      lines(x_seq, gx, col = colors[k], lwd = 3, lty = k)
    }
    
    legend("topright", legend = paste0(gs_labels, " (sigma = ", round(sigma_medians, 1), "m)"),
           col = colors, lwd = 3, lty = 1:K_classes, bg = "white", box.col = "gray80", cex = 0.9)
    
    dev.off()
    cat(sprintf("  Saved clean detection plot: %s\n", plot_png_rank))
  }
}

# Generate Multi-Species Comparative Figure (Rank 1 Models)
cat("\nGenerating combined multi-species comparison diagram plot...\n")
comb_plot_png <- file.path(detection_dir, "Detection_Plot_All_Species_Comparison.png")
png(comb_plot_png, width = 3200, height = 2400, res = 300)

par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1.5))
x_seq <- seq(0, dist_limit, length.out = 200)

for (sp_code in names(species_info)) {
  sp_name <- species_info[[sp_code]]
  r1_summary <- rank1_summary_list[[sp_code]]
  
  if (is.null(r1_summary)) next
  
  sig_rows <- r1_summary[grepl("^sigma\\[", r1_summary$Parameter), ]
  K_classes <- nrow(sig_rows)
  sig_medians <- sort(sig_rows$Median)
  
  colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")
  if (K_classes != 5) colors <- rainbow(K_classes, s = 0.8, v = 0.7)
  
  plot(NULL, xlim = c(0, dist_limit), ylim = c(0, 1.05),
       xlab = "Distance (m)", ylab = "Detection Prob g(x)",
       main = sprintf("%s (%s - Top Model)", sp_name, sp_code),
       cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.2, font.main = 2,
       las = 1, bty = "l")
  grid(col = "gray85", lty = "dotted")
  
  for (k in 1:K_classes) {
    sig <- sig_medians[k]
    gx <- exp(- (x_seq^2) / (2 * sig^2))
    lines(x_seq, gx, col = colors[k], lwd = 2.5, lty = k)
  }
  legend("topright", legend = paste0("GS ", 1:K_classes, " (", round(sig_medians, 1), "m)"),
         col = colors, lwd = 2, lty = 1:K_classes, bg = "white", cex = 0.75)
}
dev.off()
cat(sprintf("Saved combined comparison plot: %s\n", comb_plot_png))

# Export Summary Tables (Excel & CSV)
if (length(all_summary_list) > 0) {
  final_detection_excel_df <- do.call(rbind, all_summary_list)
  
  excel_path <- file.path(detection_dir, "Detection_Parameters_Summary.xlsx")
  csv_path <- file.path(detection_dir, "Detection_Parameters_Summary.csv")
  
  sheet_list <- list("All Delta2 Models Detection" = final_detection_excel_df)
  for (sp_code in names(species_info)) {
    sp_name <- species_info[[sp_code]]
    sp_df <- final_detection_excel_df %>% filter(Species_Code == sp_code)
    if (nrow(sp_df) > 0) {
      sheet_list[[sp_name]] <- sp_df
    }
  }
  
  write_xlsx(sheet_list, path = excel_path)
  write.csv(final_detection_excel_df, csv_path, row.names = FALSE)
  
  cat("\n========================================================================\n")
  cat(sprintf("SUCCESS: Exported clean detection parameter statistics for all Delta2 models:\n  - %s\n", excel_path))
  cat(sprintf("  - %s\n", csv_path))
  cat("========================================================================\n")
}
