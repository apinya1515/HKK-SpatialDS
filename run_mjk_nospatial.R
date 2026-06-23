# run_mjk_nospatial.R
library(nimble)
library(dplyr)
library(coda)
library(sf)
library(spdep)
library(parallel)

options(nimbleVerbose = FALSE)

data_tr <- read.table('line_data.txt', sep='\t', header=T)
poly <- st_read('shp/HKK1sqkmGrid.shp', quiet=TRUE)
poly$ID <- seq(1:nrow(poly))
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]
data_prop <- read.csv('TRidentity.csv')

covar_names <- colnames(data_land)[-1]
n_covar <- length(covar_names)

cat("Starting NoSpatial MCMC for Muntjac...\n")

if(!dir.exists("Results/Posteriors")) dir.create("Results/Posteriors", recursive=TRUE)

final_results <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
mjk_models <- final_results %>% filter(Species == "Muntjac" & Type == "NoSpatial")

models_to_run <- list()
for (i in 1:nrow(mjk_models)) {
  global_rank <- mjk_models$Rank[i]
  covars_str <- mjk_models$Covariates[i]
  
  rds_file <- sprintf("Results/Posteriors/Samples_MJK_NoSpatial_GlobalRank%d.rds", global_rank)
  
  w <- rep(0, n_covar)
  if (covars_str != "Intercept-only") {
    active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
    w[covar_names %in% active_vars] <- 1
  }
  models_to_run[[length(models_to_run) + 1]] <- list(mask = w, name = covars_str, rank = global_rank, file = rds_file)
}

sp_code <- "MJK"
nrep <- 10
dist_limit <- 100
data_tr_sp <- data_tr
data_tr_sp$P.dist[data_tr_sp$P.dist > dist_limit] <- dist_limit
data_sub_sp <- data_tr_sp %>% filter(Species == sp_code)
gs_max_sp <- max(data_sub_sp$Gz.sz)
gsBreaks <- unique(pmin(c(1, 2, 3, 4, 8), gs_max_sp))
dist_class_n <- 5

cl <- makeCluster(min(4, length(models_to_run)))

clusterExport(cl, c("sp_code", "nrep", "dist_limit", "gsBreaks", 
                    "dist_class_n", "data_tr_sp", "poly", "data_prop", 
                    "data_land_orig", "data_land"))
                    
clusterEvalQ(cl, {
  library(nimble)
  library(dplyr)
  library(coda)
  library(sf)
  library(spdep)
  
  species <- sp_code
  data_tr <- data_tr_sp
  
  garbage_out <- capture.output({
    source('@data_prepare_011025.R')
    source('CAR-Kumar2021-NoSpatial.R')
    
    distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
    mcmcConf <- configureMCMC(distanceModel, monitors = tracked_var, enableWAIC = TRUE)
    mcmcConf$removeSamplers('w')
    distanceMCMC <- buildMCMC(mcmcConf)
    CdistanceModel <- compileNimble(distanceModel)
    Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)
  })
  
  evaluate_model_worker <- function(model_item) {
    tryCatch({
      CdistanceModel$w <- model_item$mask
      CdistanceModel$calculate()
      
      garbage_mcmc <- capture.output({
        samples_chains <- runMCMC(Cmcmc, niter = 50000, nburnin = 30000, thin = 2, nchains = 1, WAIC=TRUE,
                                  setSeed = 999 + model_item$rank)
      })
      
      saveRDS(samples_chains, model_item$file)
      return(list(rank = model_item$rank, success = TRUE, error = NULL))
    }, error = function(e) {
      return(list(rank = model_item$rank, success = FALSE, error = e$message))
    })
  }
})

results <- parLapply(cl, models_to_run, function(m) evaluate_model_worker(m))

for (res in results) {
  if (res$success) {
    cat(sprintf("-> Successfully completed Global Rank %d\n", res$rank))
  } else {
    cat(sprintf("-> FAILED Global Rank %d: %s\n", res$rank, res$error))
  }
}

stopCluster(cl)

# After MCMC, let's extract the covariates and plot
cat("Starting MJK NoSpatial extraction and plotting...\n")
library(ggplot2)
library(tidyr)
library(ggrepel)

final_results <- read.csv("Results/Final_Model_Comparison.csv", stringsAsFactors = FALSE)
quants <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
quant_names <- paste0(quants * 100, "%")
plot_quants <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
plot_quant_names <- paste0(plot_quants * 100, "%")

for (i in 1:nrow(final_results)) {
  if (final_results$Species[i] == "Muntjac" && final_results$Type[i] == "NoSpatial") {
    rank <- final_results$Rank[i]
    covars_str <- final_results$Covariates[i]
    rds_file <- sprintf("Results/Posteriors/Samples_MJK_NoSpatial_GlobalRank%d.rds", rank)
    
    if (file.exists(rds_file)) {
      s_obj <- readRDS(rds_file)
      samps <- s_obj$samples
      
      t_quants <- quantile(samps[, "TOTAL_ABUND"], probs = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975))
      final_results$CI_2.5[i] <- t_quants["2.5%"]
      final_results$CI_5[i]   <- t_quants["5%"]
      final_results$CI_25[i]  <- t_quants["25%"]
      final_results$CI_50[i]  <- t_quants["50%"]
      final_results$CI_75[i]  <- t_quants["75%"]
      final_results$CI_95[i]  <- t_quants["95%"]
      final_results$CI_97.5[i]<- t_quants["97.5%"]
      
      active_vars <- c()
      if (covars_str != "Intercept-only") {
        active_vars <- trimws(unlist(strsplit(covars_str, "\\+")))
      }
      
      plot_df_list <- list()
      if ("beta0" %in% colnames(samps)) {
        b0_vals <- samps[, "beta0"]
        b0_q <- quantile(b0_vals, probs = quants)
        for (j in 1:length(quants)) final_results[i, paste0("beta0_", quant_names[j])] <- b0_q[j]
        plot_df_list[[length(plot_df_list) + 1]] <- data.frame(Variable = "Intercept", Value = b0_vals)
      }
      
      for (v_idx in 1:length(covar_names)) {
        v_name <- covar_names[v_idx]
        if (v_name %in% active_vars) {
          col_beta <- paste0("beta[", v_idx, "]")
          if (col_beta %in% colnames(samps)) {
            b_vals <- samps[, col_beta]
            b_q <- quantile(b_vals, probs = quants)
            for (j in 1:length(quants)) final_results[i, paste0(v_name, "_", quant_names[j])] <- b_q[j]
            plot_df_list[[length(plot_df_list) + 1]] <- data.frame(Variable = v_name, Value = b_vals)
          }
        }
      }
      
      if (length(plot_df_list) > 0) {
        plot_df <- do.call(rbind, plot_df_list)
        plot_df$Variable <- factor(plot_df$Variable, levels = c("Intercept", active_vars))
        
        quant_labels_df <- plot_df %>% group_by(Variable) %>% reframe(Quantile = plot_quant_names, Value = quantile(Value, probs = plot_quants))
        
        p <- ggplot(plot_df, aes(x = Variable, y = Value)) +
          geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.6) +
          geom_text_repel(data = quant_labels_df, aes(label = paste(Quantile, round(Value, 2), sep=": "), y = Value), size = 3, direction = "y", segment.color = "grey50", hjust = 0, nudge_x = 0.3) +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
          theme_minimal() +
          labs(title = sprintf("Muntjac - NoSpatial Rank %d", rank), subtitle = paste("Covariates:", covars_str), x = "Covariate", y = "Posterior Distribution") +
          coord_cartesian(xlim = c(1, length(unique(plot_df$Variable)) + 0.8))
        
        ggsave(sprintf("Results/Covariates/Posterior_Covariates_MJK_NoSpatial_Rank%d.jpg", rank), plot = p, width = 10, height = 7, dpi = 300)
      }
    }
  }
}
write.csv(final_results, "Results/Final_Model_Comparison.csv", row.names = FALSE)
cat("Completed.\n")
