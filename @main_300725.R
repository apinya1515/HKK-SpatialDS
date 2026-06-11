###### Based on Kumar, 2021
###### 12/04/2025 ######
### Add binary indicator ("w") w/ Bernoulli for each regression coefficient
### Add group size classes (defined by "gsBreaks" in INITIALIZATION)
### ## calculate prob of each groups size ("gs_m") and prob of each group size class ("gs_k")
### ## "gs_k" is used to calculate "pi" (product between "gs_k" and cond.prob of each distance class ("mn_cell[k, j]"))
### ## expected size of each group size class ("gs_k_mean") is used to calculate "sigma"
### Add post-training MCMC diagnostics; R-hat, ESS, and autocorrelation (for nchains >= 2 only)
###### 22/04/2025 ######
### Change custom half-normal detection prob to Nimble-compatibility function pnorm (with 2 multiplication)
### Calculate cumulative prob at largest ditance (F_dist_limit[k]) and normalized other distance class with this
### Change prior of sigma0 to lognormal to allow larger values

###### 22/07/2025
### Split code into 2 files: 
###   main_run_model = main file for pre-model and post-model
###   model*** = nimble model code to be called by source() function

###### 30/09/2025
### See model code for changes in model
### Split data preparation path into @data_prepare_[date].R

##### 5/10/2025
### Add replication in calculation of transect-level abundance
### lam_tr = nrep * (prop * lam), so, log(lam_tr) = log(nrep) + log(prop * lam)

library(nimble)
library(dplyr)
library(coda) # for post-training mcmc check 
library(MCMCvis)
library(sf) # for import spatial data
library(spdep) # for defining spatial neighbors (used for CAR)
library(terra) # for raster export and visualization
library(arm)

##############################################################################################
# Transect data (table 1)
data_tr <- read.table('line_data.txt', sep='\t', header=T)
#data_tr$Tr.no <- paste0(data_tr$Tr.no,'-',data_tr$Walk.no)

# shapefile of study area grid
poly <- st_read('shp\\HKK1sqkmGrid.shp')
poly$ID <- seq(1:nrow(poly))  # define ID directly to data instead
plot(poly)
# landscape data extracted to each grid
data_land_orig <- read.csv('HKK_Cov1sqkm_.csv')

### Covariates for model
#data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_med', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE', 'HE')]
data_land <- data_land_orig[,c('grid_id', 'dist_str', 'ndvi_cv', 'elev', 'slope', 'BB', 'DD', 'DE')]


# transect X landscape data
# Show how long each transect overlap to each grid
data_prop <- read.csv('TRidentity.csv')

###########################################################################################
########################## INITIALIZATION ####################################
### Input parameters
species <- "BTG"
# Map species code to full name for dynamic plots and filenames
species_names <- c(BTG = "Banteng", SBR = "Sambar deer", GAR = "Gaur", MJK = "Muntjac")
species_name <- if (species %in% names(species_names)) species_names[species] else species
# visualize data
discrete.hist(data_tr[data_tr$Species==species,]$Gz.sz, main='Grp size', xlab='Group Size')
hist(data_tr[data_tr$Species==species,]$P.dist, main='Distance', xlab='Distance', breaks=10)

nrep <- 10 # number of replications for each transect
dist_limit <- 100 # max observed distance
data_tr$P.dist[data_tr$P.dist > dist_limit ] <- dist_limit  # set maximum distance to the limit
gsBreaks <- c(1, 2, 3, 4, 8) # upper breaks group size of each gs classes - c(1, 3) -> 2 classes 1st class = 1, 2nd class >= 2, the maximum group size is 3
dist_class_n <- 5 # number of dist classes
##############################################################################

### Call data preparation
source('@data_prepare_011025.R')

### *****Call model file*****
## Model with indicator
#source(paste0(path, 'model-detection-function-grsize-landscape-CAR-Kumar2021-30092025.R'))
## Model without indicator
source('CAR-Kumar2021-NoSpatial.R')

####################################################################################

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, constants = constants, data = data, inits = inits)
# see initialization
distanceModel$initializeInfo()

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = tracked_var, enableWAIC = TRUE)
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

###################################################################################################
# Run the MCMC
t_start <- Sys.time() # start time
samples_chains <- runMCMC(Cmcmc, niter = 150000, nburnin = 100000, thin = 2, nchains = 3, WAIC=TRUE,
                          setSeed = c(999, 111, 555))
(t_elapse <- Sys.time() - t_start) # time elapse
###summary(samples_chains)

##############################################
# #### Continuous running
# Cmcmc$run(niter = 5000) # Run for 5000 iter first
# # Then run another more 20000 iterations
# Cmcmc$run(niter = 20000, reset = FALSE) 
# # Then run another 20000
# Cmcmc$run(niter = 20000, reset = FALSE) 

# as.list(Cmcmc$mvSamples)

################################################

### SAVE MODEL TO RDS files (So you dont have to run the same model again)
saveRDS(samples_chains, paste0(species_name, '_Model_dist_str_ndvi_cv_elev_slope_BB_DD_DE.rds'))
#samples_chains <- readRDS('Banteng_FullModel.rds')
#samples_chains <- readRDS("D:/Model_SBR_HKKcut20000burn15000.rds")
####################################################################################################



####################################################################################################
################################ Post training ####################################################
# Convergence check - use package coda
# Diagnostics
# Convert samples to coda format
mcmc.list <- as.mcmc.list(lapply(samples_chains$samples, as.mcmc))
# Filter parameters with prefixes "beta", "sigma", "w", and "muc"
param_names <- colnames(mcmc.list[[1]])
selected_params <- param_names[grep("^(beta|sigma|muc|p|sigma_spatial|sigma_alpha)", param_names)]
colnames(covar) # see covariates names
# Subset mcmc.list to include only selected parameters
sub_mcmc.list <- mcmc.list[, selected_params, drop = FALSE]

# Diagnostics plots and summaries
MCMCtrace(sub_mcmc.list, pdf = FALSE) # Visualize MCMC trace

# Wrap convergence diagnostic checks in tryCatch to prevent crashing on constant/non-varying nodes
tryCatch({
  print(gelman.diag(sub_mcmc.list))  # R-hat (R-hat < 1.1 for convergence)
}, error = function(e) {
  message("Warning: gelman.diag failed (often due to constant parameter chains like indicators): ", e$message)
})

print(effectiveSize(sub_mcmc.list))  # ESS (ESS 100 - 200 for convergence)
autocorr.plot(sub_mcmc.list)  # Autocorrelation (Low Autocorrelation for convergence)
print(MCMCsummary(sub_mcmc.list))  # Summary with R-hat and ESS

### merge multiple chains into 1
samples <- do.call(rbind, samples_chains$samples)

### Detection function parameters
if (any(grepl("^sigma\\[", colnames(samples)))) {
  sigma_all <- samples[, grep("^sigma\\[", colnames(samples)), drop = FALSE]
  sigma_k <- apply(sigma_all, 2, median)
  print("Median sigma for each group size class:")
  print(sigma_k)
  p <- median(samples[, 'p'])
  print(paste("Median p:", p))
  
  # multinomial pi
  if (any(grepl("^pi\\[", colnames(samples)))) {
    pi_val <- apply(samples[, grep("^pi\\[", colnames(samples)), drop = FALSE], 2, median)
    pi_mat <- matrix(pi_val, nrow = gs_class_n, byrow = FALSE)
    print("Multinomial detection probability matrix (pi):")
    print(pi_mat)
    print(paste("Sum of pi:", sum(pi_mat)))
    
    ### Visualize pi = joint detection prob
    plot(pi_mat[1,], type = 'n', xlab = 'Distance Class', ylab = 'Joint Probability (pi)', 
         main = 'Multinomial Detection Probability (pi) by Group Size Class',
         ylim = c(0, max(pi_mat) * 1.1))
    for(i in 1:nrow(pi_mat)){
      lines(pi_mat[i,], lty = i, col = i, lwd = 2)
    }
    legend('topright', legend = 1:nrow(pi_mat), title = 'Group size class', lty = 1:nrow(pi_mat), col = 1:nrow(pi_mat), lwd = 2)
  }
  
  # calculate and plot conditional detection prob (g(x)) given group size class
  x_seq <- seq(0, dist_limit)
  hp <- hist(data_sub$P.dist, main = 'Histogram of Observed Distances & Estimated Detection Functions', 
             xlab = 'Distance', xlim = c(0, dist_limit))
  for(i in 1:length(sigma_k)){
    gx <- exp((-1 * x_seq^2) / (2 * sigma_k[i]^2))
    lines(x_seq, gx * max(hp$counts), col = i, lwd = 2)
  }
  legend('topright', legend = 1:length(sigma_k), title = 'Group size class', col = 1:length(sigma_k), lty = 1, lwd = 2)
}

# Beta regression coefficients
beta_cols <- grep("^beta\\[", colnames(samples))
if (length(beta_cols) > 0) {
  # Sort beta columns numerically to avoid alphabetical sorting issues (e.g. beta[10] before beta[2])
  beta_names <- colnames(samples)[beta_cols]
  beta_indices <- as.numeric(gsub("beta\\[([0-9]+)\\]", "\\1", beta_names))
  beta_cols <- beta_cols[order(beta_indices)]
  beta <- samples[, beta_cols, drop = FALSE]
  
  if (length(beta_cols) == ncol(covar)) {
    colnames(beta) <- colnames(covar)
  } else {
    warning("Number of beta columns in MCMC (", length(beta_cols), ") does not match number of covariates (", ncol(covar), ")")
  }
  
  # Adjust margins to ensure vertical axis labels are not clipped
  op <- par(mar = c(7, 4, 4, 2) + 0.1)
  boxplot(beta, outline = FALSE, main = paste(species_name, "Beta Regression Coefficients"), las = 2, ylab = "Effect Size")
  abline(h = 0, lty = 3, col = "red")
  par(op) # Restore margins
  
  # Indicator variables (variable selection)
  w_cols <- grep("^w\\[", colnames(samples))
  if (length(w_cols) > 0) {
    # Sort w columns numerically
    w_names <- colnames(samples)[w_cols]
    w_indices <- as.numeric(gsub("w\\[([0-9]+)\\]", "\\1", w_names))
    w_cols <- w_cols[order(w_indices)]
    w <- samples[, w_cols, drop = FALSE]
    
    if (length(w_cols) == ncol(covar)) {
      colnames(w) <- colnames(covar)
    }
    w_prop <- colSums(w) / nrow(w)
    
    op <- par(mar = c(7, 4, 4, 2) + 0.1)
    barplot(w_prop, main = paste(species_name, "Inclusion Probability (w)"), las = 2, ylab = "Probability", ylim = c(0, 1))
    abline(h = 0.5, lty = 2, col = "gray")
    par(op) # Restore margins
    
    # Plot coefficients of variables with inclusion prob > 0.3
    active_vars <- which(w_prop > 0.3)
    if (length(active_vars) > 0) {
      beta_w <- beta[, active_vars, drop = FALSE]
      op <- par(mar = c(7, 4, 4, 2) + 0.1)
      boxplot(beta_w, outline = FALSE, main = paste(species_name, "Beta Coefficients (Inclusion > 0.3)"), las = 2, ylab = "Effect Size")
      abline(h = 0, lty = 3, col = "red")
      par(op) # Restore margins
    }
  }
}

# Grid-specific spatial random effects (b_spatial or spatial_z)
spatial_cols <- grep("^(b_spatial|spatial_z)\\[", colnames(samples))
if (length(spatial_cols) > 0) {
  b_spatial <- samples[, spatial_cols, drop = FALSE]
  bspat_quant <- apply(b_spatial, 2, function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))
  
  # Plot spatial random effect distributions
  par(mfrow = c(1, 2))
  hist(bspat_quant[3,], main = 'Median of Spatial Random Effects (CAR)', xlab = 'Random Effect')
  hist(bspat_quant[5,], main = 'Upper (95%) of Spatial Random Effects', xlab = 'Random Effect')
  par(mfrow = c(1, 1))
  
  bspat_sd <- apply(b_spatial, 2, sd, na.rm = TRUE)
  hist(bspat_sd, main = 'SD of Spatial Random Effects (CAR)', xlab = 'SD')
  
  # Spatial visualization
  poly$bspat_med <- bspat_quant[3,]
  plot(poly["bspat_med"], main = 'Median of Spatial Random Effects (CAR)', breaks = 'quantile')
}

# Grid-specific expected group abundance (z)
z_cols <- grep("^z\\[", colnames(samples))
if (length(z_cols) > 0) {
  z <- samples[, z_cols, drop = FALSE]
  z_quant <- apply(z, 2, function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))
  
  par(mfrow = c(2, 2))
  hist(z_quant[3,], main = 'Median Grid Group-Abundance', xlab = 'z', breaks = 50)
  hist(z_quant[2,], main = 'Lower (5%) Grid Group-Abundance', xlab = 'z', breaks = 50)
  hist(z_quant[4,], main = 'Upper (95%) Grid Group-Abundance', xlab = 'z', breaks = 50)
  hist(z_quant, main = 'All Grid Group-Abundance Quantiles', xlab = 'z', breaks = 100, xlim = c(0, max(z_quant[3,]) * 1.5))
  par(mfrow = c(1, 1))
  
  z_sd <- apply(z, 2, sd, na.rm = TRUE)
  hist(z_sd, main = 'SD of Grid Group-Abundance', xlab = 'z_sd', breaks = 50)
  
  poly$z_med <- z_quant[3,]
  poly$z_upper <- z_quant[4,]
  poly$z_lower <- z_quant[2,]
  poly$z_sd <- z_sd
  
  plot(poly["z_med"], main = 'Median Grid Group Abundance (z)', breaks = 'quantile')
  plot(poly["z_lower"], main = 'Lower (5%) Grid Group Abundance (z)', breaks = 'quantile')
  plot(poly["z_upper"], main = 'Upper (95%) Grid Group Abundance (z)', breaks = 'quantile')
  plot(poly["z_sd"], main = 'SD Grid Group Abundance (z)', breaks = 'quantile')
}

# Site(transect)-specific intercept (alpha)
alpha_cols <- grep("^alpha", colnames(samples))
if (length(alpha_cols) > 0) {
  alpha <- samples[, alpha_cols, drop = FALSE]
  hist(alpha, breaks = 50, main = 'Site-Specific Intercept (alpha)', xlab = 'alpha')
}

# Site-specific proportionated grid abundance (Z)
Z_cols <- grep("^Z\\[", colnames(samples))
if (length(Z_cols) > 0) {
  Z <- samples[, Z_cols, drop = FALSE]
  hist(apply(Z, 2, median), breaks = 50, main = 'Transect Proportionated Grid Abundance (Z)', xlab = 'Z')
}

# Transect-specific expected abundance (lam)
lam_cols <- grep("^lam\\[", colnames(samples))
if (length(lam_cols) > 0) {
  lam <- samples[, lam_cols, drop = FALSE]
  boxplot(lam, main = 'Transect Expected Abundance (lam)', ylab = 'lam')
}

# Average Group Size (AGS)
if ("AGS" %in% colnames(samples)) {
  avg_gs <- samples[, "AGS"]
  hist(avg_gs, main = paste('Posterior of', species_name, 'Average Group Size (AGS)'), xlab = 'Group Size', breaks = 50)
}

# Grid-specific individual Abundance (ABUND = z * AGS)
abund_cols <- grep("^ABUND\\[", colnames(samples))
if (length(abund_cols) == 0) {
  # If ABUND is not tracked, calculate it post hoc as the product of group abundance (z) and average group size (AGS)
  z_cols <- grep("^z\\[", colnames(samples))
  if (length(z_cols) > 0 && "AGS" %in% colnames(samples)) {
    z_samples <- samples[, z_cols, drop = FALSE]
    ags_samples <- samples[, "AGS"]
    abund <- sweep(z_samples, 1, ags_samples, FUN = "*")
    colnames(abund) <- paste0("ABUND[", 1:ncol(abund), "]")
    samples <- cbind(samples, abund)
    abund_cols <- grep("^ABUND\\[", colnames(samples))
  }
}

if (length(abund_cols) > 0) {
  abund <- samples[, abund_cols, drop = FALSE]
  colnames(abund) <- 1:grid_n
  
  abund_med <- apply(abund, 2, median, na.rm = TRUE)
  abund_upper <- apply(abund, 2, function(x) quantile(x, 0.95))
  abund_lower <- apply(abund, 2, function(x) quantile(x, 0.05))
  abund_sd <- apply(abund, 2, sd, na.rm = TRUE)
  
  hist(abund_med, breaks = 100, xlim = c(0, max(abund_med) * 1.5), probability = TRUE, main = 'Grid Individual Abundance Histogram', xlab = 'Abundance')
  hist(abund_upper, breaks = 100, add = TRUE, col = rgb(1,0,0,0.4), probability = TRUE)
  hist(abund_lower, breaks = 100, add = TRUE, col = rgb(0,1,0,0.4), probability = TRUE)
  legend('topright', legend = c('Median', 'Upper (95%)', 'Lower (5%)'), fill = c('white', rgb(1,0,0,0.4), rgb(0,1,0,0.4)))
  
  hist(abund_sd, main = 'SD of Grid-Level Individual Abundance', xlab = 'SD', breaks = 100)
  
  poly$abund_med <- abund_med
  poly$abund_upper <- abund_upper
  poly$abund_lower <- abund_lower
  poly$abund_sd <- abund_sd
  poly$abund_lowsd <- abund_sd < 50
  
  plot(poly["abund_med"], main = 'Median Grid-Level Individual Abundance')
  plot(poly["abund_upper"], main = 'Upper (95%) Grid-Level Individual Abundance')
  plot(poly["abund_lower"], main = 'Lower (5%) Grid-Level Individual Abundance')
  plot(poly["abund_sd"], main = 'SD Grid-Level Individual Abundance')
  plot(poly["abund_lowsd"], main = 'Abundance - Low SD Grids (SD < 100)')
  
  # Calculate total abundance post hoc from grid-wise abundance (fixing original typos)
  print("Post hoc calculated Total Abundance (sum of grid-level posteriors):")
  print(paste("Median sum:", sum(abund_med)))
  print(paste("Upper sum:", sum(abund_upper)))
  print(paste("Lower sum:", sum(abund_lower)))
  
  # Landscape density distribution
  abund_med_land <- apply(abund, 1, median, na.rm = TRUE)
  hist(abund_med_land, breaks = 100, main = 'Posterior of Median Landscape Grid Abundance', xlab = 'Median Grid Abundance')
  density_quant <- quantile(abund_med_land, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
  print("Landscape-wide Grid Abundance Quantiles:")
  print(density_quant)
  abline(v = density_quant[1], lty = 2, col = "blue")
  abline(v = density_quant[7], lty = 2, col = "blue")
  abline(v = density_quant[4], lwd = 2, col = "red")
  legend('topright', legend = c('Median', '95% CI'), col = c('red', 'blue'), lty = c(1, 2), lwd = c(2, 1))
}

# Total abundance (from within-model derived-quantities)
if ("TOTAL_ABUND" %in% colnames(samples)) {
  total_abund <- samples[, "TOTAL_ABUND"]
  
  # Calculate summary stats on the full posterior
  med_val <- median(total_abund)
  ci_vals <- quantile(total_abund, c(0.05, 0.95))
  
  # Print values to console
  ta_quant <- quantile(total_abund, c(0.05, 0.25, 0.5, 0.75, 0.95))
  print("Within-model Total Abundance Posterior Quantiles:")
  print(ta_quant)
  
  # Filter to a reasonable range for visualization if there are outliers
  plot_data <- total_abund[total_abund < 10000]
  if (length(plot_data) == 0) plot_data <- total_abund
  
  h <- hist(plot_data, breaks = 200, main = paste(species_name, 'Total Abundance Posterior'), xlab = 'Abundance')
  
  # Draw vertical lines for median and 90% CI (5% and 95% quantiles)
  abline(v = med_val, lty = 2, col = "red", lwd = 2)
  abline(v = ci_vals, col = 'blue', lty = 2, lwd = 1.5)
  
  # Add text labels for median and 90% CI bounds on the plot itself
  y_text <- max(h$counts) * 0.85
  xlims <- range(plot_data)
  
  if (med_val >= xlims[1] && med_val <= xlims[2]) {
    text(x = med_val, y = y_text, labels = sprintf("Median: %.1f", med_val), col = "red", pos = 4, cex = 0.9)
  }
  if (ci_vals[1] >= xlims[1] && ci_vals[1] <= xlims[2]) {
    text(x = ci_vals[1], y = y_text * 0.75, labels = sprintf("5%%: %.1f", ci_vals[1]), col = "blue", pos = 2, cex = 0.8)
  }
  if (ci_vals[2] >= xlims[1] && ci_vals[2] <= xlims[2]) {
    text(x = ci_vals[2], y = y_text * 0.75, labels = sprintf("95%%: %.1f", ci_vals[2]), col = "blue", pos = 4, cex = 0.8)
  }
  
  # Legend showing numeric values of Median and 90% CI
  legend_labels <- c(
    sprintf("Median: %.1f", med_val),
    sprintf("90%% CI: [%.1f, %.1f]", ci_vals[1], ci_vals[2])
  )
  legend('topright', legend = legend_labels, col = c('red', 'blue'), lty = c(2, 2), lwd = c(2, 1.5), bg = "white")
}


