# scratch/calculate_weights.R
# Post-process model selection summary files to calculate wAIC weights for the top 5 models

library(dplyr)

post_process_weights <- function(file_path) {
  if (!file.exists(file_path)) {
    cat("File not found:", file_path, "\n")
    return(NULL)
  }
  
  df <- read.csv(file_path)
  
  # Calculate wAIC weights within each species
  df_updated <- df %>%
    group_by(Species) %>%
    mutate(
      Weight = exp(-0.5 * Delta_wAIC) / sum(exp(-0.5 * Delta_wAIC))
    ) %>%
    ungroup()
  
  # Reorder columns to ensure Weight is right after Delta_wAIC
  df_updated <- df_updated[, c("Species", "Rank", "Covariates", "wAIC", "Delta_wAIC", "Weight", "lppd", "pWAIC")]
  
  write.csv(df_updated, file_path, row.names = FALSE)
  cat("Successfully updated wAIC weights in:", file_path, "\n")
}

cat("Updating model selection summary files...\n")
post_process_weights("Results/Model_Selection_Summary.csv")
post_process_weights("Results/Model_Selection_Spatial_Summary.csv")
