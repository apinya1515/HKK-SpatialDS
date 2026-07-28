# save_delta2_models.R
# Filters all models with Delta_wAIC <= 2 across species and saves:
# 1. Models_Delta2_Summary.csv & .xlsx in Results/Delta2_Models/
# 2. Copies/Saves the full MCMC .rds samples into Results/MCMC/ and Results/Delta2_Models/

library(dplyr)
library(writexl)

cat("========================================================================\n")
cat("EXTRACTING & SAVING ALL MODELS WITH Delta_wAIC <= 2\n")
cat("========================================================================\n")

output_dir <- "Results/Delta2_Models"
mcmc_dir <- "Results/MCMC"
posteriors_dir <- "Results/Posteriors"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(mcmc_dir)) dir.create(mcmc_dir, recursive = TRUE)

# Read model comparison summary
if (!file.exists("Results/Final_Model_Comparison.csv")) {
  stop("Results/Final_Model_Comparison.csv not found.")
}

final_results <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)

# Filter for models with Delta_wAIC <= 2
delta2_models <- final_results %>% 
  filter(Delta_wAIC <= 2) %>% 
  arrange(Species, Rank)

cat(sprintf("Found %d models with Delta_wAIC <= 2 across all species:\n", nrow(delta2_models)))
print(delta2_models[, c("Species", "Rank", "Type", "Covariates", "wAIC", "Delta_wAIC", "Weight")])

# Save summary tables (CSV & XLSX)
csv_out <- file.path(output_dir, "Models_Delta2_Summary.csv")
xlsx_out <- file.path(output_dir, "Models_Delta2_Summary.xlsx")

write.csv(delta2_models, csv_out, row.names = FALSE)

# Create sheet per species for Excel export
species_list <- unique(delta2_models$Species)
sheet_list <- list("All Delta2 Models" = delta2_models)

for (sp in species_list) {
  sheet_list[[sp]] <- delta2_models %>% filter(Species == sp)
}

write_xlsx(sheet_list, path = xlsx_out)
cat(sprintf("\nSaved Delta <= 2 summary to:\n  - %s\n  - %s\n", csv_out, xlsx_out))

# Map species names to species codes
all_species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

# Copy/Ensure MCMC RDS files exist in Results/MCMC/ for all Delta_wAIC <= 2 models
cat("\nChecking & copying MCMC .rds files for Delta_wAIC <= 2 models...\n")

for (i in 1:nrow(delta2_models)) {
  sp_name <- delta2_models$Species[i]
  sp_code <- all_species_codes[[sp_name]]
  rank <- delta2_models$Rank[i]
  covars <- delta2_models$Covariates[i]
  m_type <- delta2_models$Type[i]
  
  target_rds_mcmc <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank))
  target_rds_delta <- file.path(output_dir, sprintf("Samples_%s_Rank%d.rds", sp_code, rank))
  
  # Candidate source RDS files
  possible_sources <- c(
    sprintf("%s/Samples_%s_GlobalRank%d.rds", posteriors_dir, sp_code, rank),
    sprintf("%s/Samples_%s_Rank%d.rds", posteriors_dir, sp_code, rank),
    sprintf("%s/MCMC_Samples_%s.rds", mcmc_dir, sp_code)
  )
  
  found_source <- NULL
  for (src in possible_sources) {
    if (file.exists(src)) {
      found_source <- src
      break
    }
  }
  
  if (!is.null(found_source)) {
    cat(sprintf("  [%s Rank %d (%s)] Copying %s -> %s\n", sp_code, rank, m_type, found_source, target_rds_mcmc))
    file.copy(found_source, target_rds_mcmc, overwrite = TRUE)
    file.copy(found_source, target_rds_delta, overwrite = TRUE)
  } else {
    cat(sprintf("  [%s Rank %d (%s)] MCMC RDS not present locally yet (will be generated upon run).\n", sp_code, rank, m_type))
  }
}

cat("\n========================================================================\n")
cat("SUCCESS: All models with Delta_wAIC <= 2 processed and saved!\n")
cat("========================================================================\n")
