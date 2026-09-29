# summarize_rerun.R
# Post-run records for the converged re-runs (run after the new RDS files are in Results/MCMC):
#   1. Results/model_summary/Rerun_MCMC_Settings.csv
#      iterations per chain, draws kept, final max Rhat / min ESS per model (for the Methods text)
#   2. Results/model_summary/Rerun_Old_vs_New_Estimates.csv
#      TOTAL_ABUND, beta0, active betas, sigma_spatial: old (archived) vs new posterior summaries
# Usage: Rscript summarize_rerun.R <archive_dir>

suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
archive_dir <- args[1]
covar_names <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')
species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")

final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
models <- final_df %>% filter(Delta_wAIC < 2) %>% arrange(Species, Rank)

settings <- list(); est <- list()
for (r in seq_len(nrow(models))) {
  m <- models[r, ]
  sp <- species_codes[[m$Species]]
  tag <- sprintf("%s_Rank%d", sp, m$Rank)

  lg <- read.csv(file.path("Results/MCMC_rerun", tag, "convergence_log.csv"), stringsAsFactors = FALSE)
  last <- lg[nrow(lg), ]
  # burn-in = first floor(k/2) of k chunks (logged directly by the current rerun_models.R)
  burnin <- if (!is.null(last$burnin_iterations) && !is.na(last$burnin_iterations)) last$burnin_iterations else
    floor(last$chunk / 2) * last$iterations_per_chain / last$chunk
  settings[[r]] <- data.frame(
    Species = m$Species, Species_Code = sp, Rank = m$Rank, Type = m$Type, Covariates = m$Covariates,
    Chains = NA_integer_, Iterations_per_chain = last$iterations_per_chain,
    Burnin_per_chain = burnin, Thin = 10,
    Draws_kept_per_chain = last$kept_draws_per_chain,
    Rule = "Rhat < 1.045 for every parameter; 30,000-100,000 iterations per chain",
    Max_Rhat_nonCAR = round(last$max_rhat_noncar, 4), Max_Rhat_CAR_nodes = round(last$max_rhat_car, 4),
    Rhat_sigma_spatial = round(last$sigma_spatial_rhat, 4),
    Min_ESS_key = round(last$min_ess_key, 0), Min_ESS_param = last$min_ess_param,
    Converged = last$pass, stringsAsFactors = FALSE
  )

  summ <- function(f) {
    s_list <- readRDS(f)
    n_chains_seen <<- length(s_list)
    s <- do.call(rbind, s_list)
    pars <- c("TOTAL_ABUND", "beta0",
              sprintf("beta[%d]", which(covar_names %in% trimws(unlist(strsplit(m$Covariates, "\\+"))))),
              "sigma_spatial")
    pars <- intersect(pars, colnames(s))
    out <- t(sapply(pars, function(p) c(Mean = mean(s[, p]), Median = median(s[, p]),
                                         q2.5 = unname(quantile(s[, p], 0.025)), q97.5 = unname(quantile(s[, p], 0.975)))))
    data.frame(Parameter = pars, out, row.names = NULL)
  }
  new <- summ(sprintf("Results/MCMC/MCMC_Samples_%s.rds", tag))
  settings[[r]]$Chains <- n_chains_seen
  old_f <- file.path(archive_dir, sprintf("MCMC_Samples_%s.rds", tag))
  old <- if (file.exists(old_f)) summ(old_f) else new[0, ]
  names(old)[-1] <- paste0("Old_", names(old)[-1]); names(new)[-1] <- paste0("New_", names(new)[-1])
  e <- merge(new, old, by = "Parameter", all.x = TRUE)
  e$Parameter <- ifelse(grepl("^beta\\[", e$Parameter),
                        paste0(e$Parameter, " ", covar_names[as.integer(gsub("\\D", "", e$Parameter))]), e$Parameter)
  e$Pct_change_median <- round(100 * (e$New_Median - e$Old_Median) / abs(e$Old_Median), 1)
  est[[r]] <- cbind(Species = m$Species, Rank = m$Rank, Type = m$Type, Covariates = m$Covariates, e)
  invisible(gc())
}
write.csv(do.call(rbind, settings), "Results/model_summary/Rerun_MCMC_Settings.csv", row.names = FALSE)
write.csv(do.call(rbind, est), "Results/model_summary/Rerun_Old_vs_New_Estimates.csv", row.names = FALSE)
cat("Saved Rerun_MCMC_Settings.csv and Rerun_Old_vs_New_Estimates.csv in Results/model_summary/\n")
