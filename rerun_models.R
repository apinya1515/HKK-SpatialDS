# rerun_models.R
# Re-runs one model of Results/Final_Model_Comparison.csv (species + rank), one chain per
# process, until the rule below is met. run_rerun_queue.ps1 runs all Delta_wAIC < 2 models.
#
# Usage (one process per chain, all N_CHAINS (default 3) chains of a model at the same time):
#   Rscript rerun_models.R <SPECIES_CODE> <RANK> <CHAIN 1..N>
# Check-only mode (no sampling): apply the rule to chunks on disk (capped at max_iter); write
# the final RDS if it passes, or with FORCE_SAVE=1 at the cap. Set CHUNK_ITER to their size.
#   Rscript rerun_models.R <SPECIES_CODE> <RANK> check
#
# Spatial: CAR-Kumar2021-IND-NC.R (non-centred, sparse; same posterior as CAR-Kumar2021-IND.R).
# NoSpatial: CAR-Kumar2021-NoSpatial.R with the rank's covariate mask applied.
#
# After each chunk of `chunk_iter`, chain 1 checks the post-burn-in draws (burn-in = first
# floor(k/2) of k chunks) and writes CONTINUE / STOP:
#   STOP when >= 30,000 iterations per chain AND Rhat < 1.045 for EVERY parameter
#   (non-spatial parameters, sigma_spatial, and all 3013 b_spatial grid nodes),
#   or at 100,000 iterations per chain (then saved and logged as NOT converged).
# ESS is logged only. Rhat = same formula as generate_model_summaries.R (xlsx).
# On STOP chain 1 writes Results/MCMC_rerun/MCMC_Samples_<SP>_Rank<N>.rds in the same format as
# Results/MCMC/*.rds: list(chain1, chain2, chain3) of matrices, plus attributes
# "iteration" (iteration number of each kept draw) and "thin".

suppressMessages({
  library(nimble); library(dplyr); library(coda); library(sf); library(spdep)
})
options(nimbleVerbose = FALSE)

args <- commandArgs(trailingOnly = TRUE)
sp_code <- args[1]
rank <- as.integer(args[2])
check_only <- args[3] == "check"
chain <- if (check_only) 1L else as.integer(args[3])

n_chains <- as.integer(Sys.getenv("N_CHAINS", "3"))
chunk_iter <- as.integer(Sys.getenv("CHUNK_ITER", "10000"))
thin <- 10
min_iter <- as.integer(Sys.getenv("MIN_ITER", "30000"))   # env overrides are for smoke tests only
max_iter <- as.integer(Sys.getenv("MAX_ITER", "100000"))
max_chunks <- max_iter %/% chunk_iter
rhat_max <- 1.045
keep_per_chain <- 5000 # draws per chain kept for diagnostics and the final file

out_dir <- Sys.getenv("RERUN_DIR", "Results/MCMC_rerun")
work_dir <- file.path(out_dir, sprintf("%s_Rank%d", sp_code, rank))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  who <- if (check_only) "check" else sprintf("chain%d", chain)
  msg <- sprintf("[%s %s R%d %s] %s", format(Sys.time(), "%m-%d %H:%M:%S"), sp_code, rank, who, sprintf(...))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(work_dir, "progress.log"), append = TRUE)
}

species_codes <- c("Banteng" = "BTG", "Sambar deer" = "SBR", "Gaur" = "GAR", "Muntjac" = "MJK", "Wild boar" = "PIG")
final_df <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
top <- final_df[final_df$Species == names(species_codes)[species_codes == sp_code] & final_df$Rank == rank, ]
stopifnot(nrow(top) == 1)
spatial <- top$Type == "Spatial"

# ---- Data (same preparation as run_compiled_posteriors.R) ----
data_tr <- read.table('line_data.txt', sep = '\t', header = TRUE)
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet = TRUE)
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[, c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')
covar_names <- colnames(data_land)[-1]

species <- sp_code
nrep <- 10
dist_limit <- 100
dist_class_n <- 5
data_tr$P.dist[data_tr$P.dist > dist_limit] <- dist_limit
gs_max_sp <- max(data_tr$Gz.sz[data_tr$Species == species], na.rm = TRUE)
gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
invisible(capture.output(source('@data_prepare_011025.R')))

active_vars <- trimws(unlist(strsplit(top$Covariates, "\\+")))
w_mask <- as.numeric(covar_names %in% active_vars)
active_beta <- sprintf("beta[%d]", which(w_mask == 1))

# ---- Model definition (derive_abundance comes from the NC model file) ----
set.seed(7919 * chain + 101 * rank + match(sp_code, species_codes))
if (spatial) {
  source('CAR-Kumar2021-IND-NC.R')
  monitors <- tracked_var
} else {
  source('CAR-Kumar2021-NoSpatial.R')
  inits$w <- w_mask
  inits$beta <- rnorm(ncol(covar), 0, 0.5)
  monitors <- c("beta0", "beta", "muc", "sigma0", "sigma", "p", "pi", "gs_k", "AGS", "TOTAL_ABUND")
}

# ---- Diagnostics helpers ----
calc_rhat_vec <- function(chains) { # same formula as generate_model_summaries.R (m chains)
  n <- nrow(chains[[1]])
  w <- Reduce(`+`, lapply(chains, function(ch) apply(ch, 2, var))) / length(chains)
  b_over_n <- apply(sapply(chains, colMeans), 1, var)
  var_hat <- ((n - 1) / n) * w + b_over_n
  ifelse(w == 0, 1.0, sqrt(pmax(1.0, var_hat / w)))
}

chunk_file <- function(k, ch) file.path(work_dir, sprintf("chunk%02d_chain%d.rds", k, ch))
decision_file <- function(k) file.path(work_dir, sprintf("decision%02d.txt", k))
wait_for <- function(f) while (!file.exists(f)) Sys.sleep(5)

# Load post-burn-in chunks of all chains, add derived columns (final format)
load_kept <- function(k) {
  first <- floor(k / 2) + 1
  draws_per_chunk <- chunk_iter / thin
  iter_all <- (first - 1) * chunk_iter + thin * seq_len((k - first + 1) * draws_per_chunk)
  sel <- unique(round(seq(1, length(iter_all), length.out = min(keep_per_chain, length(iter_all)))))
  ch <- lapply(seq_len(n_chains), function(chn) {
    m <- do.call(rbind, lapply(first:k, function(j) readRDS(chunk_file(j, chn))))
    stopifnot(nrow(m) == length(iter_all))
    m <- m[sel, , drop = FALSE]
    if (spatial) {
      d <- derive_abundance(m, covar, w_mask, water_mask)
      colnames(d$b_spatial) <- sprintf("b_spatial[%d]", seq_len(ncol(d$b_spatial)))
      m <- cbind(m[, !grepl("^u\\[", colnames(m)), drop = FALSE], TOTAL_ABUND = d$TOTAL_ABUND, d$b_spatial)
    }
    m[, sort(colnames(m))]
  })
  list(chains = ch, iteration = iter_all[sel], burnin = (first - 1) * chunk_iter)
}

check_convergence <- function(k) {
  kept <- load_kept(k)
  ch <- kept$chains
  pn <- colnames(ch[[1]])
  is_car <- grepl("^b_spatial\\[", pn)
  rh <- calc_rhat_vec(ch)
  rh[!is.finite(rh)] <- Inf # non-finite = not converged
  key <- intersect(c("beta0", active_beta, "sigma0", "p", "muc", "AGS", "TOTAL_ABUND", "sigma_spatial"), pn)
  ess <- sapply(key, function(p) tryCatch(unname(coda::effectiveSize(coda::mcmc.list(lapply(ch, function(m) coda::mcmc(m[, p]))))), error = function(e) NA_real_))
  res <- data.frame(
    chunk = k, chunk_iter = chunk_iter, iterations_per_chain = k * chunk_iter, burnin_iterations = kept$burnin,
    kept_draws_per_chain = nrow(ch[[1]]),
    max_rhat_noncar = max(rh[!is_car]), worst_noncar = pn[!is_car][which.max(rh[!is_car])],
    max_rhat_car = if (any(is_car)) max(rh[is_car]) else NA,
    sigma_spatial_rhat = if ("sigma_spatial" %in% pn) rh[pn == "sigma_spatial"] else NA,
    min_ess_key = suppressWarnings(min(ess, na.rm = TRUE)), min_ess_param = c(names(ess)[which.min(ess)], NA)[1],
    stringsAsFactors = FALSE
  )
  res$pass <- k * chunk_iter >= min_iter && max(rh) < rhat_max
  list(res = res, kept = kept)
}

log_check <- function(r) { # append (tolerates logs written by earlier versions of this script)
  f <- file.path(work_dir, "convergence_log.csv")
  old <- if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
  write.csv(bind_rows(old, r), f, row.names = FALSE)
  log_msg("check: max Rhat non-CAR %.4f (%s) | max Rhat CAR %s | sigma_spatial %s | min ESS %.0f (%s) -> %s",
          r$max_rhat_noncar, r$worst_noncar, format(round(r$max_rhat_car, 4)), format(round(r$sigma_spatial_rhat, 4)),
          r$min_ess_key, r$min_ess_param, ifelse(r$pass, "CONVERGED", "continue"))
}

save_final <- function(chk, converged) {
  out <- setNames(chk$kept$chains, paste0("chain", seq_len(n_chains)))
  attr(out, "iteration") <- chk$kept$iteration
  attr(out, "thin") <- thin
  saveRDS(out, file.path(out_dir, sprintf("MCMC_Samples_%s_Rank%d.rds", sp_code, rank)))
  log_msg("saved final samples (%s)", ifelse(converged, "converged", "MAX CHUNKS REACHED - NOT converged"))
}

# ---- Check-only mode ----
if (check_only) {
  k <- 0
  while (all(file.exists(chunk_file(k + 1, seq_len(n_chains))))) k <- k + 1
  if (k == 0) stop("no complete chunk sets in ", work_dir)
  k <- min(k, max_chunks) # never beyond the iteration cap
  force <- Sys.getenv("FORCE_SAVE") == "1" && k == max_chunks
  chk <- check_convergence(k)
  log_check(chk$res)
  if (chk$res$pass || force) save_final(chk, chk$res$pass)
  cat(ifelse(chk$res$pass, "PASS", "FAIL"), "\n")
  quit(save = "no", status = 0)
}

# ---- Build and compile ----
model <- nimbleModel(distanceModelCode, constants = constants, data = data, inits = inits, check = FALSE)
if (!spatial) model$w <- w_mask
model$calculate()

conf <- configureMCMC(model, monitors = monitors, thin = thin)
if (!spatial) conf$removeSamplers('w')
# Extra RW_block on intercept + active slopes (on top of the defaults)
conf$addSampler(target = c("beta0", active_beta), type = "RW_block")
if (spatial) {
  conf$removeSamplers("sigma_spatial")
  conf$addSampler(target = "sigma_spatial", type = "slice")
}
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(buildMCMC(conf), project = model)
log_msg("Compiled %s model: %s (w = %s)", top$Type, top$Covariates, paste(w_mask, collapse = ""))

#   
# ---- Chunked sampling ----
for (k in 1:max_chunks) {
  Cmcmc$run(chunk_iter, reset = (k == 1), resetMV = TRUE, progressBar = FALSE)
  smp <- as.matrix(Cmcmc$mvSamples)
  tmp <- paste0(chunk_file(k, chain), ".tmp")
  saveRDS(smp, tmp)
  file.rename(tmp, chunk_file(k, chain))
  rm(smp); invisible(gc())
  log_msg("chunk %d done (%d iterations)", k, k * chunk_iter)

  if (chain == 1) {
    for (other in seq_len(n_chains)[-1]) wait_for(chunk_file(k, other))
    chk <- check_convergence(k)
    log_check(chk$res)
    stop_now <- chk$res$pass || k == max_chunks
    if (stop_now) save_final(chk, chk$res$pass)
    rm(chk); invisible(gc())
    writeLines(ifelse(stop_now, "STOP", "CONTINUE"), decision_file(k))
  } else {
    wait_for(decision_file(k))
    stop_now <- readLines(decision_file(k))[1] == "STOP"
  }
  if (stop_now) break
}
log_msg("finished")
