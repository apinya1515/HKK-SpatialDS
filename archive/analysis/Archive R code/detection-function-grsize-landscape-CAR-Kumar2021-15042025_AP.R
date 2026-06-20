#### Based on Kumar, 2021
### 12/04/2025
### Add binary indicator ("ind") w/ Bernoulli for each regression coefficient
### Add interval group size classes (defined by "gsBreaks" in INITIALIZATION)
### ## calculate prob of each groups size ("gs_m") and prob of each class ("gs_k")
### ## "gs_k" is used to calculate "pi" (product between "gs_k" and cond.prob of each distance class ("mn_cell[k, j]"))
### ## mean of each group size class ("gs_k_mean") is used to calculate "sigma"

library(nimble)
library(dplyr)
setwd("D:/2024SpatialDistance/analysis")
# Transect data (table 1)
#data_tr <- read.table('D:\\UserData\\Dropbox\\@KU_pop_course\\day1_exercise\\line_data.txt', 
#                       sep='\t', header=T)

data_tr <- read.csv('D:\\2024SpatialDistance\\analysis\\hkkdata.csv')#Apinya
tran_name <- sort(unique(data_tr$Tr.no)) # name of transect
ntran <- length(tran_name) # number of transect

# landscape data
#data_land <- read.csv('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\HKK_Cov1sqkm_.csv')
data_land <- read.csv('D:\\2024SpatialDistance\\analysis\\HKK_Cov1sqkm_15042025.csv')#Apinya
grid_n <- nrow(data_land) # number of grid
grid_id <- data_land$grid_id # grid id

# select only covariate rows and rescale
covar <- as.data.frame(data_land[,-1] %>% scale())

# transect X landscape data
#data_prop <- read.csv('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\TRidentity.csv')
data_prop <- read.csv('D:\\2024SpatialDistance\\analysis\\TRidentity.csv')#Apinya
##############################################################################
########################## INITIALIZATION ####################################
### Input parameters
species <- 'SBR'
dist_limit <- 120 # max observed distance
dist_class_n <- 2 # number of dist classes
gsBreaks <- c(1, 3) # upper breaks group size of each gs classes - c(1, 3) = 2 classes 1st = 1, 2nd >= 2
##############################################################################

### Create observation matrix (matrix=transect, row=cluster size class ,column=distance class)
dist_width <- dist_limit/dist_class_n # interval between distance class
distBreaks = seq(dist_width, dist_class_n*dist_width, by=dist_width) # upper breaks of distance classes

# subset data by species
data_sub <- data_tr %>% filter(Species==species)

gs_max <- max(data_sub$Gz.sz)   # max group size
gs_class_n <- length(gsBreaks) # number og gs classes

#### Create zeros-array for the data
y_matrix <- array(0, dim=c(gs_class_n, dist_class_n, ntran), 
                  dimnames=list(1:gs_class_n, distBreaks, tran_name)) # 3 dimensions = group size, distance, transect

### Populate y data for target species ... by +1 to y_matrix for each data point that fall into the cell
# y_matrix[group_size, distance, transect]
for(i in 1:nrow(data_sub)){
  samp <- data_sub[i,] # get data
  if(samp$P.dist > dist_limit){print(paste0('data ',i, ' distance:',  round(samp$P.dist,2), ' OUT OF RANGE !!!'));next;}
  tr <- which(tran_name ==  samp$Tr.no) # transect index
  j <- sum(distBreaks <= samp$P.dist) + 1 # distance class index
  k <- min(which(gsBreaks >= samp$Gz.sz)) # group size class index
  print(paste0('data ',i, ' gs:', samp$Gz.sz, ' distance:',  round(samp$P.dist,2),' in grsz:', k, ' dist:', j, ' transect:', tr))
  y_matrix[k, j, tr] = y_matrix[k, j, tr] + 1 # add observation into group class - dist class - transect
}
print(y_matrix)

sum(y_matrix) # check the number of observations
nrow(data_sub) # check the number of observations

### Grid - Transect Mapping ###
# Create a matrix row=transect col=grid to map proportion from each grid to transect
head(data_prop)
lenM <- matrix(0, nrow=ntran, ncol=grid_n) # create zeo matrix row=#of trasect col=#of grid
for(i in 1:nrow(data_prop)){
  dat <- data_prop[i,]
  lenM[dat$tr_seq, dat$grid_id] <- dat$Length
}
# total length for each transect
lenSum <- rowSums(lenM)
# create proportion matrix
propM <- sweep(lenM, 1, lenSum, FUN="/")
propM[1:30, 1:30]
rowSums(propM) # check that each row sum to 1

### Spatial dependency ######################
### CAR elements #############################
library(sf)
library(spdep)
#poly <- st_read('D:\\UserData\\KUDrive\\@Projects\\HKK_Line_Transect\\shp\\HKK1sqkmGrid.shp')
poly <- st_read('D:\\2024SpatialDistance\\analysis\\HKK1sqkmGrid.shp')#Apinya
poly$ID <- seq(1:grid_n)
# Create a neighbors list (Queen's contiguity)
nb <- poly2nb(poly, queen = TRUE, row.names = poly$ID)
adj <- unlist(nb) # adjacency index (grid id of each neighbor)
njoin <- length(adj) # number of spatial joins
num <- sapply(nb, length) # number of neighbors of each grid 
### (sum of num must match with the length of adj)

####################################################################################
####################################################################################
####################################################################################
# Model code in NIMBLE
distanceModelCode <- nimbleCode({
  ##### Priors #####
  for(i in 1:n_covar){
    beta[i] ~ dnorm(0, 0.001) # priors of landscape var regression coefficient
    ind[i] ~ dbern(0.5) # priors of binary indicator for each coefficient
  }
  muc ~ dunif(1, 3) # average cluster size
  sigma0 ~ dunif(0.1, 5) # sigma0
  p ~ dunif(0.01, 10) # model parameter for sigma ~ group_size
  
  # CAR priors
  tau ~ dgamma(0.001, 0.001) # precision param for spatial CAR
  weights[1:njoin] <- 1 # weight of all spatial joins = 1
  b_spatial[1:L] ~ dcar_normal(adj = adj[1:njoin], weights = weights[1:njoin], 
                           num = num[1:L], tau = tau, zero_mean = 0) # random spatial effect (intercept) intercept for each grid
  
  ##### Model #####
  
  ### Landscape abundance for each grid
  ind_beta[1:n_covar] <- ind[1:n_covar] * beta[1:n_covar]
  for(l in 1:L) {
    # Regression abundance + b_spatial(CAR elements)
    log_lam[l] <- inprod(covar[l, 1:n_covar], ind_beta[1:n_covar]) + b_spatial[l] # regression for log abundance (without intercept)
    lam[l] <- exp(log_lam[l]) # lambda for each grid
  }
  
  ### Abundance at each transect
  for (i in 1:I) {
    tr_lam[i] <- inprod(propM[i,1:L], lam[1:L])
  }
  
  ### Model prob group size class
  for(m in 1:gs_max) { # loop through every group size
    ### calculate truncated Poisson (at m=1) prob for each size class from 1 to gs_max
    gs_m[m] <- (exp(-muc) * pow(muc, m)) / (factorial(m) * (1 - exp(-muc)))
  }
  ### calculate prob of each group size class & mean group size for each class
  if(gsBreaks[1] == 1) { # if the first class contain only group size (m) = 1
    gs_k[1] <- gs_m[1]  # Direct assignment for single element
    gs_k_mean[1] <- gs_m[1] * 1 / gs_k[1]  # Mean for class 1
  } else {
    gs_k[1] <- sum(gs_m[1:gsBreaks[1]])  # Sum for multiple elements
    gs_k_mean[1] <- sum(gs_m[1:gsBreaks[1]] * 1:gsBreaks[1]) / gs_k[1]  # Mean for class 1 (see page 45)
  }
  for(k in 2:K) { # loop through group size class 2 and on
    gs_k[k] <- sum(gs_m[(gsBreaks[k-1]+1):gsBreaks[k]]) # sum of gs_m for class 2 and on
    gs_k_mean[k] <- sum(gs_m[(gsBreaks[k-1]+1):gsBreaks[k]] * (gsBreaks[k-1]+1):gsBreaks[k]) /  gs_k[k] # gs mean for class 2 and on (see page 45)
  }
  
  ### Calculate sigma for each size class
  # sigma is function of group size class (k) bigger group size -> larger sigma[k] -> slower decay half normal detection function
  for(k in 1:K){
    log(sigma[k]) <- sigma0 + p * (gs_k_mean[k] - 1)
  }
  
  ### Distance sampling model
  # Compute pi for each distance class [j] and group size class [k]
  # pi[k.j] = probability of detection is depend on group size(k) and distance(j) but is the same across transect
  for(k in 1:K){ # loop for distance class = 1
    # mn_cell = conditional multinomial prob for distance detection function j of group size class k
    # mn_cell[k,1] is for halfnormal at first distance class
    # (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) = scaling of standard normal fn
    # phi(distBreaks[1] = cumulative prob at class j=1
    # 0.5 = the minus part to remove cumulative prob of the left half of standard normal fn
    mn_cell[k, 1] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * 
      (phi(distBreaks[1] / sigma[k]) - 0.5)
    # p[k,j] = probability that the obs is in class k and j
    pi[k, 1] <- mn_cell[k, 1] * gs_k[k]
    for(j in 2:J){ # loop for distance class j >= 2
      # do the same but remove the cum prob of previous dist class j-1 instead
      mn_cell[k, j] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * 
        (phi(distBreaks[j] / sigma[k]) - phi(distBreaks[j - 1] / sigma[k]))
      # p[k,j]
      pi[k, j] <- mn_cell[k, j] * gs_k[k]
    }
  }

  #### Model data for each transect
  # Loop for each element in all data matrix
  for(i in 1:I){
    for(k in 1:K){
      for(j in 1:J){
        # expected observation in transect i, group size class k, and distance class j
        mu[k, j, i] <- tr_lam[i] * pi[k, j] # probability of class k and j * abundance of transect i
        ### Model transect data
        y[k, j, i] ~ dpois(mu[k, j, i])
      }
    }
  }
  
})

####################################################################################

constants <- list(I = ntran, J = dist_class_n, 
                  K = gs_class_n, gs_max = gs_max,
                  L = grid_n, n_covar = ncol(covar),
                  distBreaks = distBreaks, gsBreaks = gsBreaks,
                  adj = adj, num = num, njoin = njoin)
data <- list(y = y_matrix, covar = covar, propM = propM)
inits <- list(muc = runif(1, 1, 3), sigma0 = runif(1, 1, 2), p = runif(1, 0.5, 1), tau = runif(1,0.5,1),
              beta = rnorm(ncol(covar), 0,1), b_spatial = runif(grid_n, 0, 0.1), 
              ind = rep(1, ncol(covar)))

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, 
                             constants = constants, 
                             data = data, 
                             inits = inits)

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = c("sigma", "p", "gs_m","gs_k", "sigma0", "pi", 'lam', 'beta', 'tau', 'b_spatial', 'ind'))
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

# Run the MCMC
t_start <- Sys.time() # start time
samples <- runMCMC(Cmcmc, niter = 1500, nburnin = 1000, thin = 5, nchains = 1)
(t_elapse <- Sys.time() - t_start) # time elapse

#summary(samples)

################################ Post training ############################
# fitted detection function
sigma_1 <- median(samples[,'sigma[1]'])
sigma_2 <- median(samples[,'sigma[2]'])
p <- median(samples[,'p'])
p

# multinomial pi
pi_val <- apply(samples[,grep("^pi", colnames(samples))], 2, median)
pi <- matrix(pi_val, nrow=gs_class_n, byrow = F)
pi
sum(pi)

### Visualize pi
plot(pi[1,], type='l', xlab='Distance Class', ylab='Group Size', main='Multinomial Prob')
lines(pi[2,], lty=2)
legend('topright', legend=c('1','2'), title = 'Group size', lty=c(1,2))

# calculate gx (only distance)
x_seq <- seq(0, dist_limit)
gx_1 <- exp((-1*x_seq^2)/(2*sigma_1^2))
gx_2 <- exp((-1*x_seq^2)/(2*sigma_2^2))
hp <- hist(data_sub$P.dist, main='Histogram & Detection function', xlab='Distance')
lines(x_seq, gx_1*max(hp$counts), col='red', lwd=2)
lines(x_seq, gx_2*max(hp$counts), col='green', lwd=2)
legend('topright', legend=c('1','2'), title = 'Group class', col=c('red', 'green'), lty=1)

# Site-specific Lambda (Abundance)
lam <- samples[,grep("^lam", colnames(samples))]
tr_names <- dimnames(y_matrix)[[3]]
colnames(lam) <- 1:grid_n
boxplot(lam[,1:50], ylim=c(0,100))
boxplot(lam[,101:150], ylim=c(0,100))
abun_med <- apply(lam, 2, median, na.rm = TRUE)
hist(abun_med)

# Site-specific random effect
b_spatial <- samples[,grep("^b_spatial", colnames(samples))]
bspat_med <- apply(b_spatial, 2, median, na.rm = TRUE)

# Beta
beta <- samples[,grep("^beta", colnames(samples))]
colnames(beta) <- colnames(covar)
boxplot(beta)
abline(h=0, lty=3)

# Indicator
ind <- samples[,grep("^ind", colnames(samples))]
colnames(ind) <- colnames(covar)
barplot(ind)

# Inject back to shp file
library(terra)
abun_med[abun_med>20]
poly$abun_med <- abun_med
poly$bspat_med <- bspat_med
plot(poly["abun_med"], main='Median grid-level abundance')
plot(poly["bspat_med"], main='Median CAR effects')

