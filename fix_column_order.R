# fix_column_order.R
# One-off repair (29/09/2026): rerun_models.R saved columns in alphabetical order
# (b_spatial[1], b_spatial[10], b_spatial[100], ...), which scrambled grid order for scripts that
# take b_spatial columns positionally (step4_model_averaging.R). Reorders every sample matrix to
# natural NIMBLE order (by variable, then numeric index), keeping all attributes.

natural_order <- function(cn) {
  base <- sub("\\[.*$", "", cn)
  inside <- ifelse(grepl("[", cn, fixed = TRUE), sub("^.*?\\[", "", cn, perl = TRUE), "")
  idx <- regmatches(inside, gregexpr("[0-9]+", inside))
  i1 <- sapply(idx, function(v) if (length(v) >= 1) as.numeric(v[1]) else 0)
  i2 <- sapply(idx, function(v) if (length(v) >= 2) as.numeric(v[2]) else 0)
  order(base, i1, i2)
}

files <- c(list.files("Results/MCMC", "^MCMC_Samples_.*\\.rds$", full.names = TRUE),
           list.files("Results/MCMC_rerun", "^MCMC_Samples_.*\\.rds$", full.names = TRUE))
for (f in files) {
  s <- readRDS(f)
  att <- attributes(s)
  o <- natural_order(colnames(s[[1]]))
  if (identical(o, seq_along(o))) { cat("already ordered:", f, "\n"); next }
  s2 <- lapply(s, function(m) m[, o, drop = FALSE])
  attributes(s2) <- att
  stopifnot(identical(colnames(s2[[1]])[grep("^b_spatial", colnames(s2[[1]]))][1:3],
                      c("b_spatial[1]", "b_spatial[2]", "b_spatial[3]")) || !any(grepl("^b_spatial", colnames(s2[[1]]))))
  saveRDS(s2, f)
  cat("reordered:", f, "\n")
}
