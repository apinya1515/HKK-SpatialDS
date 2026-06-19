# model_selection.R
# Backward stepwise variable selection using wAIC in NIMBLE
# For each species, parallelized model evaluations across 8 cores

library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)
library(parallel)

# Set global options
options(nimbleVerbose = FALSE)

# Load data in main process
data_tr <- read.table('line_data.txt', sep='\t', header=T)
poly <- st_read('shp/HKK1sqkmGrid.shp')
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')

# Define species and covariates
all_species <- c("BTG", "SBR", "GAR", "MJK")
species_names <- c(BTG = "Banteng", SBR = "Sambar deer", GAR = "Gaur", MJK = "Muntjac")
covar_names <- colnames(data_land)[-1] # dist_str, ndvi_cv, elev, slope, BB, DD, DE
n_covar <- length(covar_names)

cat("Starting model variable selection based on wAIC using backward stepwise search...\n")
cat("Number of covariates: ", n_covar, " (", paste(covar_names, collapse=", "), ")\n\n")

# Store final results for all species
all_results <- list()

for (sp in all_species) {
  sp_name <- species_names[sp]
  cat("========================================================================\n")
  cat("SPECIES: ", sp_name, " (", sp, ")\n")
  cat("========================================================================\n")
  
  # Configure species parameters
  species <- sp
  nrep <- 10
  dist_limit <- 100
  
  # Filter data for this species to get max group size
  data_sub_sp <- data_tr %>% filter(Species == sp)
  data_tr_sp <- data_tr
  data_tr_sp$P.dist[data_tr_sp$P.dist > dist_limit] <- dist_limit
  gs_max_sp <- max(data_sub_sp$Gz.sz)
  
  # Breaks (dynamically adjusted to prevent out-of-bounds errors for small group size species)
  gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
  dist_class_n <- 5
  
  cat("Initializing cluster of 4 workers for parallel model runs to prevent OOM...\n")
  cl <- makeCluster(4)
  
  # Export variables to workers
  clusterExport(cl, c("sp", "species_names", "nrep", "dist_limit", "gsBreaks", 
                      "dist_class_n", "data_tr_sp", "poly", "data_prop", 
                      "data_land_orig", "data_land"))
  
  # Initialize environment on workers
  clusterEvalQ(cl, {
    library(nimble)
    library(dplyr)
    library(coda)
    library(sf)
    library(spdep)
    
    # Set species specific variables
    species <- sp
    data_tr <- data_tr_sp
    
    # Prepare data for this worker (suppress excessive console output to prevent parallel deadlock)
    garbage_out <- capture.output({
      source('@data_prepare_011025.R')
      
      # Load model
      source('CAR-Kumar2021-NoSpatial.R')
      
      # Build model
      distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
      
      # Configure MCMC
      mcmcConf <- configureMCMC(distanceModel, monitors = tracked_var, enableWAIC = TRUE)
      mcmcConf$removeSamplers('w')
      
      # Build and compile
      distanceMCMC <- buildMCMC(mcmcConf)
      CdistanceModel <- compileNimble(distanceModel)
      Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)
    })
    
    # Worker function for running one model
    evaluate_model_worker <- function(w_mask) {
      tryCatch({
        CdistanceModel$w <- w_mask
        CdistanceModel$calculate() # Propagate w values
        
        # Suppress MCMC progress output to prevent deadlock
        garbage_mcmc <- capture.output({
          samples_chains <- runMCMC(Cmcmc, niter = 50000, nburnin = 30000, thin = 2, nchains = 1, WAIC=TRUE,
                                    setSeed = 999)
        })
        
        list(
          w = w_mask,
          WAIC = samples_chains$WAIC$WAIC,
          lppd = samples_chains$WAIC$lppd,
          pWAIC = samples_chains$WAIC$pWAIC,
          error = NULL
        )
      }, error = function(e) {
        list(
          w = w_mask,
          WAIC = Inf,
          lppd = NA,
          pWAIC = NA,
          error = e$message
        )
      })
    }
  })
  
  cat("Compilation completed on all workers.\n\n")
  
  # Stepwise tracking variables
  evaluated_models <- list()
  
  # Helper to print and log evaluated model
  log_model <- function(w_mask, waic, lppd, pwaic, error = NULL) {
    covars_str <- if (sum(w_mask) == 0) "Intercept-only" else paste(covar_names[w_mask == 1], collapse = " + ")
    w_str <- paste(w_mask, collapse = "")
    
    # Check if already evaluated
    if (w_str %in% names(evaluated_models)) {
      return(evaluated_models[[w_str]])
    }
    
    model_res <- list(
      covariates = covars_str,
      w = w_mask,
      WAIC = waic,
      lppd = lppd,
      pWAIC = pwaic,
      error = error
    )
    evaluated_models[[w_str]] <<- model_res
    
    if (is.null(error)) {
      cat(sprintf("Model: %-50s | wAIC: %8.2f | pWAIC: %6.2f\n", covars_str, waic, pwaic))
    } else {
      cat(sprintf("Model: %-50s | Error: %s\n", covars_str, error))
    }
    return(model_res)
  }
  
  # 1. Backward stepwise selection starting from the Full Model
  cat("--- STARTING BACKWARD STEPWISE SELECTION ---\n")
  cat("Running Full Model...\n")
  w_curr <- rep(1, n_covar)
  res_full <- parLapply(cl, list(w_curr), function(w) evaluate_model_worker(w))[[1]]
  curr_res_back <- log_model(w_curr, res_full$WAIC, res_full$lppd, res_full$pWAIC, res_full$error)
  
  w_curr_back <- w_curr
  step <- 1
  converged <- FALSE
  
  while (!converged) {
    active_vars <- which(w_curr_back == 1)
    if (length(active_vars) == 0) {
      cat("No covariates left to drop.\n")
      break
    }
    
    cat(sprintf("\nBackward Step %d: Evaluating models with one variable dropped (current best wAIC: %.2f)...\n", step, curr_res_back$WAIC))
    
    candidate_masks <- list()
    uncached_masks <- list()
    for (v in active_vars) {
      w_cand <- w_curr_back
      w_cand[v] <- 0
      w_str <- paste(w_cand, collapse = "")
      candidate_masks[[length(candidate_masks) + 1]] <- w_cand
      if (!(w_str %in% names(evaluated_models))) {
        uncached_masks[[length(uncached_masks) + 1]] <- w_cand
      }
    }
    
    if (length(uncached_masks) > 0) {
      results_cand <- parLapply(cl, uncached_masks, function(w) evaluate_model_worker(w))
      for (res in results_cand) {
        log_model(res$w, res$WAIC, res$lppd, res$pWAIC, res$error)
      }
    }
    
    logged_results <- list()
    for (w_cand in candidate_masks) {
      w_str <- paste(w_cand, collapse = "")
      logged_results[[length(logged_results) + 1]] <- evaluated_models[[w_str]]
    }
    
    best_cand_idx <- which.min(sapply(logged_results, function(x) x$WAIC))
    best_cand_res <- logged_results[[best_cand_idx]]
    
    if (best_cand_res$WAIC < curr_res_back$WAIC) {
      cat(sprintf("-> Dropping '%s' improved wAIC from %.2f to %.2f.\n", 
                  covar_names[w_curr_back == 1 & best_cand_res$w == 0], curr_res_back$WAIC, best_cand_res$WAIC))
      w_curr_back <- best_cand_res$w
      curr_res_back <- best_cand_res
      step <- step + 1
    } else {
      cat("-> No single variable drop improved wAIC. Backward search finished.\n")
      converged <- TRUE
    }
  }
  
  # 2. Forward stepwise selection starting from Intercept-only Model
  cat("\n--- STARTING FORWARD STEPWISE SELECTION ---\n")
  cat("Running Intercept-only Model...\n")
  w_curr <- rep(0, n_covar)
  w_str <- paste(w_curr, collapse = "")
  # Check if intercept-only model is already run
  if (w_str %in% names(evaluated_models)) {
    curr_res_for <- evaluated_models[[w_str]]
  } else {
    res_null <- parLapply(cl, list(w_curr), function(w) evaluate_model_worker(w))[[1]]
    curr_res_for <- log_model(w_curr, res_null$WAIC, res_null$lppd, res_null$pWAIC, res_null$error)
  }
  
  w_curr_for <- w_curr
  step <- 1
  converged <- FALSE
  
  while (!converged) {
    inactive_vars <- which(w_curr_for == 0)
    if (length(inactive_vars) == 0) {
      cat("All covariates are already active.\n")
      break
    }
    
    cat(sprintf("\nForward Step %d: Evaluating models with one variable added (current best wAIC: %.2f)...\n", step, curr_res_for$WAIC))
    
    candidate_masks <- list()
    uncached_masks <- list()
    for (v in inactive_vars) {
      w_cand <- w_curr_for
      w_cand[v] <- 1
      w_str <- paste(w_cand, collapse = "")
      candidate_masks[[length(candidate_masks) + 1]] <- w_cand
      if (!(w_str %in% names(evaluated_models))) {
        uncached_masks[[length(uncached_masks) + 1]] <- w_cand
      }
    }
    
    if (length(uncached_masks) > 0) {
      results_cand <- parLapply(cl, uncached_masks, function(w) evaluate_model_worker(w))
      for (res in results_cand) {
        log_model(res$w, res$WAIC, res$lppd, res$pWAIC, res$error)
      }
    }
    
    logged_results <- list()
    for (w_cand in candidate_masks) {
      w_str <- paste(w_cand, collapse = "")
      logged_results[[length(logged_results) + 1]] <- evaluated_models[[w_str]]
    }
    
    best_cand_idx <- which.min(sapply(logged_results, function(x) x$WAIC))
    best_cand_res <- logged_results[[best_cand_idx]]
    
    if (best_cand_res$WAIC < curr_res_for$WAIC) {
      cat(sprintf("-> Adding '%s' improved wAIC from %.2f to %.2f.\n", 
                  covar_names[w_curr_for == 0 & best_cand_res$w == 1], curr_res_for$WAIC, best_cand_res$WAIC))
      w_curr_for <- best_cand_res$w
      curr_res_for <- best_cand_res
      step <- step + 1
    } else {
      cat("-> No single variable addition improved wAIC. Forward search finished.\n")
      converged <- TRUE
    }
  }
  
  # 3. Choose the lowest wAIC model
  cat("\n--- SEARCH SUMMARY ---\n")
  cat(sprintf("Backward selection model wAIC: %.2f (%s)\n", curr_res_back$WAIC, curr_res_back$covariates))
  cat(sprintf("Forward selection model wAIC: %.2f (%s)\n", curr_res_for$WAIC, curr_res_for$covariates))
  
  if (curr_res_back$WAIC <= curr_res_for$WAIC) {
    cat(sprintf("-> Selected model from Backward Stepwise (wAIC: %.2f).\n\n", curr_res_back$WAIC))
  } else {
    cat(sprintf("-> Selected model from Forward Stepwise (wAIC: %.2f).\n\n", curr_res_for$WAIC))
  }
  
  cat("Shutting down parallel cluster...\n")
  stopCluster(cl)
  
  # Process and sort all evaluated models for this species
  df_models <- do.call(rbind, lapply(evaluated_models, function(x) {
    data.frame(
      Covariates = x$covariates,
      wAIC = x$WAIC,
      lppd = x$lppd,
      pWAIC = x$pWAIC,
      stringsAsFactors = FALSE
    )
  }))
  
  # Remove failed models
  df_models <- df_models %>% filter(is.finite(wAIC))
  
  # Sort by wAIC
  df_models <- df_models %>% arrange(wAIC)
  
  # Calculate delta wAIC
  best_waic <- df_models$wAIC[1]
  df_models$Delta_wAIC <- df_models$wAIC - best_waic
  
  # Calculate wAIC model weights (Akaike weights equivalent for wAIC)
  df_models$Weight <- exp(-0.5 * df_models$Delta_wAIC) / sum(exp(-0.5 * df_models$Delta_wAIC))
  
  # Keep top 5 models OR any model with Delta_wAIC <= 2
  top_5_df <- df_models %>% filter(row_number() <= 5 | Delta_wAIC <= 2)
  top_5_df$Rank <- 1:nrow(top_5_df)
  top_5_df$Species <- sp_name
  
  # Reorder columns
  top_5_df <- top_5_df[, c("Species", "Rank", "Covariates", "wAIC", "Delta_wAIC", "Weight", "lppd", "pWAIC")]
  
  all_results[[sp]] <- top_5_df
  
  cat(sprintf("\nTop 5 models for %s:\n", sp_name))
  print(top_5_df)
  cat("\n\n")
}

cat("========================================================================\n")
cat("SUMMARY OF ALL SPECIES MODEL SELECTION\n")
cat("========================================================================\n")

# Combine results into one data frame
final_summary <- do.call(rbind, all_results)
print(final_summary)

# Write results to CSV file
write.csv(final_summary, "Results/Model_Selection_Summary.csv", row.names = FALSE)
cat("\nResults successfully saved to 'Results/Model_Selection_Summary.csv'\n")
