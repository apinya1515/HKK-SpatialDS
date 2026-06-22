# update_plots_retroactive.R
library(nimble)

cat("Starting retroactive plot updating...\n")

all_species_codes <- c("BTG" = "Banteng", "SBR" = "Sambar deer", "GAR" = "Gaur", "MJK" = "Muntjac")

# Get list of all Samples RDS files
rds_files <- list.files("Results/Posteriors", pattern = "^Samples_.*\\.rds$", full.names = TRUE)

for (f in rds_files) {
  fname <- basename(f)
  cat("Processing:", fname, "\n")
  
  # Parse filename to extract sp_code and suffix (e.g. Rank1 or GlobalRank1)
  # Expected format: Samples_<sp_code>_<suffix>.rds
  match_res <- regmatches(fname, regexec("^Samples_([A-Z]+)_(GlobalRank[0-9]+|Rank[0-9]+)\\.rds$", fname))[[1]]
  if (length(match_res) < 3) {
    cat("  Skipping (pattern mismatch)\n")
    next
  }
  
  sp_code <- match_res[2]
  suffix <- match_res[3]
  sp_name <- all_species_codes[[sp_code]]
  
  # Format suffix for title (e.g. "GlobalRank1" -> "Global Rank 1", "Rank1" -> "Rank 1")
  suffix_title <- suffix
  suffix_title <- gsub("GlobalRank", "Global Rank ", suffix_title)
  suffix_title <- gsub("Rank", "Rank ", suffix_title)
  
  # Load samples
  s_obj <- readRDS(f)
  samps <- s_obj$samples
  total_abund_samples <- samps[, "TOTAL_ABUND"]
  
  # Calculate CI
  quants <- quantile(total_abund_samples, probs = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975))
  
  # Define plot output path
  plot_jpg <- sprintf("Results/Posteriors/Posterior_%s_%s.jpg", sp_code, suffix)
  
  # Generate corrected plot
  jpeg(plot_jpg, width=800, height=600)
  
  x_max <- as.numeric(quants["97.5%"]) * 1.15
  plot_data <- total_abund_samples[total_abund_samples <= x_max]
  
  hist(plot_data, breaks=50, 
       main=sprintf("Total Abundance Posterior: %s (%s)", sp_name, suffix_title), 
       xlab="Total Abundance (Truncated at 97.5th %ile * 1.15)", col="lightgray", border="white", 
       xlim=c(min(plot_data), x_max))
  
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
  cat("  Saved corrected plot:", plot_jpg, "\n")
  
  # Clean up memory
  rm(s_obj, samps, total_abund_samples, plot_data)
  gc()
}

cat("Finished updating all individual posterior plots!\n")
