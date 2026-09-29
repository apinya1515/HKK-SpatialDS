# compute_waic_posthoc.R
# WAIC for each Delta_wAIC < 2 model, computed from the saved MCMC samples in Results/MCMC.
# Matches NIMBLE's default (conditional) WAIC: every data node y[k, j, i] is its own observation,
#   lppd  = sum_n log(mean_s p(y_n | theta_s))
#   pWAIC = sum_n var_s(log p(y_n | theta_s))
#   WAIC  = -2 * (lppd - pWAIC)
# mu[k, j, i] = lam[i] * pi[k, j] is rebuilt from beta0, beta (masked by the model's covariates),
# b_spatial at the transect grids and pi, exactly as in the model code.
#
# Output: Results/tables/WAIC_Rerun_Comparison.csv (old wAIC from Final_Model_Comparison.csv
# next to the recomputed one). Final_Model_Comparison.csv itself is NOT modified.

suppressMessages({ library(dplyr); library(sf); library(spdep) })

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")
covar_names <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')

# Pointwise log-likelihood matrix (draws x data nodes) for one model's samples
pointwise_loglik <- function(samps, sp_code, covars_str) {
  data_tr <- read.table('line_data.txt', sep = '\t', header = TRUE)
  poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE); poly$ID <- seq(1:nrow(poly))
  data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
  data_land <- data_land_orig[, c('grid_id', covar_names)]
  data_prop <- read.csv('TRidentity.csv')
  species <- sp_code; nrep <- 10; dist_limit <- 100; dist_class_n <- 5
  data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
  gs_max_sp <- max(data_tr$Gz.sz[data_tr$Species == species], na.rm = TRUE)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  invisible(capture.output(source('@data_prepare_011025.R', local = TRUE)))

  w_mask <- as.numeric(covar_names %in% trimws(unlist(strsplit(covars_str, "\\+"))))
  gT <- which(colSums(propM) > 0)
  covarT <- as.matrix(covar)[gT, , drop = FALSE]
  area_tr <- (2 * distBreaks[dist_class_n] * lenSum) / 1e6

  beta <- samps[, sprintf("beta[%d]", 1:7), drop = FALSE]
  lin <- tcrossprod(sweep(beta, 2, w_mask, "*"), covarT) + samps[, "beta0"]   # draws x T
  b_cols <- sprintf("b_spatial[%d]", gT)
  if (all(b_cols %in% colnames(samps))) lin <- lin + samps[, b_cols, drop = FALSE]
  zT <- sweep(exp(lin), 2, water_mask[gT], "*")
  lam <- sweep(tcrossprod(zT, propM[, gT, drop = FALSE]), 2, nrep * area_tr, "*")  # draws x I

  K <- gs_class_n; J <- dist_class_n; I <- tran_n
  ll <- matrix(NA_real_, nrow(samps), K * J * I)
  col <- 0
  for (i in 1:I) for (j in 1:J) for (k in 1:K) {
    col <- col + 1
    mu <- lam[, i] * samps[, sprintf("pi[%d, %d]", k, j)]
    ll[, col] <- dpois(y_matrix[k, j, i], mu, log = TRUE)
  }
  ll
}

waic_from_loglik <- function(ll) {
  m <- apply(ll, 2, max)
  lppd <- sum(m + log(colMeans(exp(sweep(ll, 2, m, "-")))))
  pwaic <- sum(apply(ll, 2, var))
  c(WAIC = -2 * (lppd - pwaic), lppd = lppd, pWAIC = pwaic)
}

if (sys.nframe() == 0) {
  final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
  models <- final_df %>% filter(Delta_wAIC < 2) %>% arrange(Species, Rank)
  rows <- list()
  for (r in seq_len(nrow(models))) {
    m <- models[r, ]
    sp_code <- species_codes[[m$Species]]
    f <- sprintf("Results/MCMC/MCMC_Samples_%s_Rank%d.rds", sp_code, m$Rank)
    s <- readRDS(f)
    samps <- do.call(rbind, s)
    w <- waic_from_loglik(pointwise_loglik(samps, sp_code, m$Covariates))
    cat(sprintf("%s Rank %d (%s): WAIC %.2f (old %.2f)\n", sp_code, m$Rank, m$Type, w["WAIC"], m$wAIC))
    rows[[r]] <- data.frame(Species = m$Species, Species_Code = sp_code, Old_Rank = m$Rank, Type = m$Type,
                            Covariates = m$Covariates, Old_wAIC = m$wAIC, Old_Delta_wAIC = m$Delta_wAIC,
                            New_wAIC = unname(w["WAIC"]), New_lppd = unname(w["lppd"]), New_pWAIC = unname(w["pWAIC"]))
    rm(s, samps); invisible(gc())
  }
  out <- do.call(rbind, rows) %>%
    group_by(Species) %>%
    mutate(New_Delta_within_rerun = New_wAIC - min(New_wAIC),
           New_Rank_within_rerun = rank(New_wAIC, ties.method = "first")) %>%
    ungroup() %>% arrange(Species, New_wAIC)
  dir.create("Results/tables", showWarnings = FALSE)
  write.csv(out, "Results/tables/WAIC_Rerun_Comparison.csv", row.names = FALSE)
  cat("Saved Results/tables/WAIC_Rerun_Comparison.csv\n")
}
