# generate_traceplots.R
# MCMC traceplots (all chains in each RDS: 3 in the re-runs) for every Delta_wAIC <= 2 model.
# Each panel = trace (left) + per-chain posterior density (right), titled with
# Rhat and ESS computed exactly as in Model_Convergence_and_Parameters_Summary.xlsx.
#
# Parameters traced per model:
#   beta0, ACTIVE beta[] (named by covariate), sigma0, p, muc, AGS, log10(TOTAL_ABUND)
#   Spatial models add: sigma_spatial, and the single b_spatial node (of 3013) with
#   the worst Rhat. (mean(b_spatial) is not traced: NIMBLE's dcar_normal sampler
#   centres the CAR effects each iteration, so it is ~1e-16 noise.)
#
# Inactive betas (w = 0) are not traced: with w fixed at 0 they have no likelihood
# contribution and are pure prior draws.
#
# Ranks whose RDS files are byte-identical (e.g. the MJK NoSpatial run copied to
# every NoSpatial rank by run_compiled_posteriors.R) get ONE plot, named for all ranks.
#
# Traces show the LAST `window_iter` (10,000) post-burn-in iterations, with the x-axis in
# MCMC iterations (from the RDS "iteration" attribute written by rerun_models.R). Rhat / ESS
# in the titles use ALL saved post-burn-in draws, so they match the xlsx.
#
# Outputs saved to Results/traceplots/:
#   Traceplot_[SP]_Rank[N].png  (or Rank[a]-[b] for shared files)
#   Traceplot_Diagnostics.csv   (Rhat / ESS of every traced parameter)

library(coda)
library(dplyr)

cat("========================================================================\n")
cat("GENERATING MCMC TRACEPLOTS FOR ALL DELTA_wAIC <= 2 MODELS\n")
cat("========================================================================\n")

trace_dir <- "Results/traceplots"
mcmc_dir <- "Results/MCMC"
window_iter <- 10000 # iterations shown in each trace
rhat_flag_at <- 1.045 # convergence threshold used by rerun_models.R
if (!dir.exists(trace_dir)) dir.create(trace_dir, recursive = TRUE)
unlink(list.files(trace_dir, pattern = "^Traceplot_.*\\.png$", full.names = TRUE)) # no stale plots

if (!file.exists("Results/Final_Model_Comparison.csv")) {
  stop("Results/Final_Model_Comparison.csv not found.")
}
final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
delta2_df <- final_df %>% filter(Delta_wAIC <= 2) %>% arrange(Species, Rank)

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")
covar_names <- c('dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE') # order of beta[1:7]

# Chain colours: categorical slots 1-3 of the validated reference palette (blue, orange, aqua)
chain_cols <- c("#2a78d6", "#eb6834", "#1baf7a")
ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
grid_col <- "#e4e3df"

# Vectorized (non-split) Gelman-Rubin Rhat for m chains (list of draws x params matrices) --
# same formula as generate_model_summaries.R
calc_rhat_vec <- function(chains) {
  n <- nrow(chains[[1]])
  w <- Reduce(`+`, lapply(chains, function(ch) apply(ch, 2, var))) / length(chains)
  means <- sapply(chains, colMeans)
  if (is.null(dim(means))) means <- matrix(means, nrow = 1)
  b_over_n <- apply(means, 1, var)
  var_hat <- ((n - 1) / n) * w + b_over_n
  ifelse(w == 0, 1.0, sqrt(pmax(1.0, var_hat / w)))
}

# Rhat shown on the plots uses the SAME formula as Model_Convergence_and_Parameters_Summary.xlsx
# (calc_rhat_vec), so every plotted value matches the xlsx. coda's gelman.diag PSRF, which adds
# the (1 + 1/m) between-chain factor and a df correction and so reads higher, goes to the CSV only.
calc_rhat_ess <- function(xs) { # xs = list of per-chain vectors
  ml <- coda::mcmc.list(lapply(xs, coda::mcmc))
  rhat <- calc_rhat_vec(lapply(xs, matrix))
  rhat_coda <- tryCatch(coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf[1, 1],
                        error = function(e) NA_real_)
  ess <- tryCatch(as.numeric(coda::effectiveSize(ml)), error = function(e) NA_real_)
  c(unname(rhat)[1], unname(ess)[1], median(unlist(xs), na.rm = TRUE), unname(rhat_coda)[1])
}

# One panel: trace on the left, sideways per-chain density on the right (shared y-axis)
draw_panel <- function(xs, it, label, rhat, ess) {
  yl <- range(unlist(xs), finite = TRUE)

  # Trace (x = MCMC iteration)
  par(mar = c(2.2, 3.6, 2.0, 0.2))
  plot(NULL, xlim = range(it), ylim = yl, xlab = "", ylab = "", axes = FALSE)
  abline(h = pretty(yl), col = grid_col, lwd = 0.8)
  for (j in seq_along(xs)) lines(it, xs[[j]], col = adjustcolor(chain_cols[j], 0.5), lwd = 0.5)
  axis(1, col = ink_muted, col.axis = ink_secondary, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.3, 0))
  axis(2, col = ink_muted, col.axis = ink_secondary, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0), las = 1)
  rhat_flag <- if (!is.na(rhat) && rhat >= rhat_flag_at) "  ▲" else ""
  title(main = sprintf("%s    Rhat = %.3f%s   ESS = %.0f", label, rhat, rhat_flag, ess),
        adj = 0, cex.main = 0.85, font.main = 2, col.main = ink_primary, line = 0.6)

  # Density (rotated)
  par(mar = c(2.2, 0.2, 2.0, 0.6))
  ds <- lapply(xs, function(x) density(x[is.finite(x)]))
  plot(NULL, xlim = c(0, max(sapply(ds, function(d) max(d$y))) * 1.05), ylim = yl, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yl), col = grid_col, lwd = 0.8)
  for (j in seq_along(ds)) lines(ds[[j]]$y, ds[[j]]$x, col = chain_cols[j], lwd = 1.5)
}

# Group ranks that point at byte-identical RDS files
delta2_df$rds_path <- file.path(mcmc_dir, sprintf("MCMC_Samples_%s_Rank%d.rds",
                                                  species_codes[delta2_df$Species], delta2_df$Rank))
delta2_df <- delta2_df[file.exists(delta2_df$rds_path), ]
delta2_df$md5 <- unname(tools::md5sum(delta2_df$rds_path))
groups <- split(seq_len(nrow(delta2_df)), paste(delta2_df$Species, delta2_df$md5))
groups <- groups[order(sapply(groups, function(g) paste(delta2_df$Species[g[1]], sprintf("%03d", delta2_df$Rank[g[1]]))))]

diag_rows <- list()

for (g in groups) {
  row <- delta2_df[g[1], ]
  sp_name <- row$Species
  sp_code <- species_codes[[sp_name]]
  ranks <- delta2_df$Rank[g]
  m_type <- row$Type
  shared <- length(g) > 1
  rank_tag <- if (shared) sprintf("Rank%d-%d", min(ranks), max(ranks)) else sprintf("Rank%d", row$Rank)

  cat(sprintf("\nProcessing [%s %s (%s)]: %s\n", sp_code, rank_tag, m_type, row$rds_path))
  s_obj <- readRDS(row$rds_path)
  chains <- unname(as.list(s_obj)) # 2 or 3 chains
  n_ch <- length(chains)
  # Iteration number of each saved draw; files without the attribute are indexed by draw
  iter <- attr(s_obj, "iteration")
  iter_known <- !is.null(iter)
  if (!iter_known) iter <- seq_len(nrow(chains[[1]]))
  win <- if (iter_known) which(iter > max(iter) - window_iter) else seq_along(iter)
  rm(s_obj); gc()
  per_chain <- function(col) lapply(chains, function(ch) ch[, col])

  # Active covariates. A shared NoSpatial file was run with w = 1 for all covariates
  # (the mask is never applied in run_compiled_posteriors.R), so every beta is active.
  if (shared) {
    active <- covar_names
  } else {
    active <- trimws(unlist(strsplit(row$Covariates, "\\+")))
  }
  beta_idx <- which(covar_names %in% active)

  # Each trace: vals = per-chain vectors to plot; raw = per-chain vectors for diagnostics
  traces <- list()
  add_trace <- function(label, col) traces[[label]] <<- list(vals = per_chain(col), raw = NULL)
  add_trace("beta0 (intercept)", "beta0")
  for (b in beta_idx) add_trace(sprintf("beta[%d]  %s", b, covar_names[b]), sprintf("beta[%d]", b))
  for (pn in c("sigma0", "p", "muc", "AGS")) add_trace(pn, pn)
  # Plotted on log10 scale, but diagnostics use raw TOTAL_ABUND to match the xlsx
  traces[["TOTAL_ABUND (log10 axis)"]] <- list(vals = lapply(per_chain("TOTAL_ABUND"), log10),
                                               raw = per_chain("TOTAL_ABUND"))

  sp_cols <- grep("^b_spatial\\[", colnames(chains[[1]]), value = TRUE)
  if (length(sp_cols) > 0) {
    add_trace("sigma_spatial", "sigma_spatial")
    sp_rhat <- calc_rhat_vec(lapply(chains, function(ch) ch[, sp_cols]))
    worst <- sp_cols[which.max(sp_rhat)]
    add_trace(sprintf("%s  (worst-Rhat grid)", worst), worst)
  }
  rm(chains); gc()

  # Diagnostics for every traced parameter
  diag <- t(sapply(traces, function(tr) calc_rhat_ess(if (!is.null(tr$raw)) tr$raw else tr$vals)))
  colnames(diag) <- c("Rhat", "ESS", "Median", "Rhat_coda")

  # Layout: 2 parameter columns, each = trace (3) + density (1); header row on top
  n_par <- length(traces)
  n_row <- ceiling(n_par / 2)
  mat <- matrix(0, nrow = n_row, ncol = 4)
  for (k in seq_len(n_par)) {
    r <- (k - 1) %/% 2 + 1
    cpos <- ((k - 1) %% 2) * 2 + 1
    mat[r, cpos:(cpos + 1)] <- c(2 * k, 2 * k + 1)
  }
  mat <- rbind(rep(1, 4), mat)

  out_png <- file.path(trace_dir, sprintf("Traceplot_%s_%s.png", sp_code, rank_tag))
  png(out_png, width = 3600, height = 300 + 620 * n_row, res = 300, bg = "#fcfcfb")
  layout(mat, widths = c(3, 1, 3, 1), heights = c(300 / 620, rep(1, n_row)))

  # Header
  par(mar = c(0, 1, 0, 1))
  plot.new()
  rank_txt <- if (shared) sprintf("Ranks %s", paste(ranks, collapse = ", ")) else sprintf("Rank %d", row$Rank)
  text(0, 0.78, sprintf("%s (%s) - %s - %s model", sp_name, sp_code, rank_txt, m_type),
       adj = 0, cex = 1.25, font = 2, col = ink_primary)
  sub_txt <- if (shared) {
    "Shared posterior: one NoSpatial run with all 7 covariates, copied to every NoSpatial rank"
  } else {
    sprintf("Covariates: %s  |  Delta_wAIC = %.2f", row$Covariates, row$Delta_wAIC)
  }
  win_txt <- if (iter_known) {
    sprintf("Iterations %s-%s shown (last %s)  |  Rhat/ESS from all %d saved draws/chain",
            format(min(iter[win]), big.mark = ","), format(max(iter), big.mark = ","),
            format(window_iter, big.mark = ","), length(iter))
  } else {
    sprintf("%d draws/chain (iteration numbers not stored)", length(iter))
  }
  text(0, 0.38, sprintf("%s  |  %s", sub_txt, win_txt), adj = 0, cex = 0.85, col = ink_secondary)
  legend(x = 0.78, y = 0.95, legend = paste("Chain", seq_len(n_ch)), col = chain_cols[seq_len(n_ch)], lwd = 2,
         bty = "n", cex = 0.85, text.col = ink_primary, horiz = TRUE, xjust = 0)
  text(1, 0.12, sprintf("▲ = Rhat ≥ %.2f", rhat_flag_at), adj = 1, cex = 0.7, col = ink_muted)

  for (k in seq_len(n_par)) {
    draw_panel(lapply(traces[[k]]$vals, function(v) v[win]), iter[win], names(traces)[k],
               diag[k, "Rhat"], diag[k, "ESS"])
  }
  dev.off()
  cat(sprintf("  Saved: %s\n", out_png))

  diag_rows[[out_png]] <- data.frame(
    Species = sp_name, Species_Code = sp_code, Ranks = paste(ranks, collapse = ";"),
    Type = m_type, Covariates = if (shared) "all 7 (shared NoSpatial run)" else row$Covariates,
    Parameter = names(traces), Median = signif(diag[, "Median"], 4),
    Rhat = round(diag[, "Rhat"], 4), ESS = round(diag[, "ESS"], 1),
    Rhat_coda_PSRF = round(diag[, "Rhat_coda"], 4),
    Plot_File = basename(out_png), stringsAsFactors = FALSE, row.names = NULL
  )
  rm(traces); gc()
}

diag_df <- do.call(rbind, diag_rows)
rownames(diag_df) <- NULL
diag_csv <- file.path(trace_dir, "Traceplot_Diagnostics.csv")
write.csv(diag_df, diag_csv, row.names = FALSE)
cat(sprintf("\nSaved diagnostics table: %s\n", diag_csv))

cat("\n========================================================================\n")
cat("TRACEPLOT GENERATION COMPLETE!\n")
cat("========================================================================\n")
