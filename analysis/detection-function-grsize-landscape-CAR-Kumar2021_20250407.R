#### This code consider group size as categorical

library(nimble)
library(dplyr)

# Transect data (table 1)
#data_tr <- read.table('D:\\UserData\\Dropbox\\@KU_pop_course\\day1_exercise\\line_data.txt', 
#                      sep='\t', header=T)
data_tr <- read.csv('D:\\2024SpatialDistance\\analysis\\hkkdata.csv')#Apinya
head(data_tr)
# landscape data
#data_land <- read.csv('D:\\UserData\\KUDrive\\@Projects\\Kluay_Line_Transect\\HKK_Cov1sqkm_.csv')
data_land <- read.csv('D:\\2024SpatialDistance\\analysis\\hkkcov1sqkm.csv')#Apinya
head(data_land)
grid_n <- nrow(data_land)
grid_id <- data_land$grid_id

covar <- as.data.frame(data_land[,-1] %>% scale())

# transect - landscape data
#data_prop <- read.csv('D:\\UserData\\KUDrive\\@Projects\\Kluay_Line_Transect\\TRidentity.csv')
data_prop <- read.csv('D:\\2024SpatialDistance\\analysis\\TRidentity.csv')#Apinya
head(data_prop)
##############################################################################
### Input parameters
species <- 'BTG'
dist_limit <- 100 # max observed distance, 90for SBR, 100 for BTG
dist_class_n <- 6 # number of dist class
##############################################################################

tran_name <- sort(unique(data_tr$Tr.no)) # name of transect
ntran <- length(tran_name) # number of transect

### Create observation matrix (matrix=transect, row=cluster size class ,column=distance class)
dist_width <- dist_limit/dist_class_n # interval between distance class
dist_class_upper <- seq(1, dist_class_n, by=1) * dist_width # upper bound cutoff between classes
dist_class_lower <- dist_class_upper - dist_width
dist_class_mid <- (dist_class_upper + dist_class_lower)/2

# subset data to identify cluster size classes
data_sub <- data_tr %>% filter(Species==species) # Sambardeer=SBR#Banteng=BTG

group_max <- max(data_sub$Gz.sz)   # max group size
group_class_size <- seq(1,group_max)  # group size in each class (1 increment)#
group_class_n <- length(group_class_size) # number of group size classes
##Apinya try to group group size by 2 for BTG 
group_class_size <- seq(1,group_max,by=2)  # group size in each class (1 increment)#

distBreaks = seq(dist_width, dist_class_n*dist_width, by=dist_width)

#### Create zeros-array for the data
y_matrix <- array(0, dim=c(group_class_n, dist_class_n, ntran), 
                  dimnames=list(group_class_size, dist_class_upper, tran_name)) # 3 dimensions = group size, distance, transect

### Populate y data for target species ... by +1 to y_matrix for each data point that fall into the cell
# y_matrix[group_size, distance, transect]
for(i in 1:nrow(data_sub)){
  samp <- data_sub[i,] # get data
  if(samp$P.dist > dist_limit){print(paste0('data ',i, ' distance:',  round(samp$P.dist,2), ' OUT OF RANGE !!!'));next;}
  tr <- which(tran_name ==  samp$Tr.no) # transect index
  j <- sum(dist_class_lower <= samp$P.dist) # distance class index
  k <- which(group_class_size == samp$Gz.sz) # group size class index
  print(paste0('data ',i, ' distance:',  round(samp$P.dist,2),' in grsz:', j, ' dist:', k, ' transect:', tr))
  y_matrix[k, j, tr] = y_matrix[k, j, tr] + 1 # add observation into group class - dist class - transect
}
print(y_matrix)

sum(y_matrix) # check the number of observations
nrow(data_sub) # check the number of observations

### Grid - Transect Mapping ###
# Create big matrix row=transect col=grid to map proportion from each grid to transect
head(data_prop)
lenM <- matrix(0, nrow=ntran, ncol=grid_n)
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
### CAR model #############################
library(sf)
library(spdep)
#poly <- st_read('D:\\UserData\\KUDrive\\@Projects\\Kluay_Line_Transect\\shp\\HKK1sqkmGrid.shp')
poly <- st_read('D:\\2024SpatialDistance\\analysis\\HKK1sqkmGrid.shp')#Apinya
poly$ID <- seq(1:grid_n)
# Create a neighbors list (Queen's contiguity)
nb <- poly2nb(poly, queen = TRUE, row.names = poly$ID)
adj <- unlist(nb) # adjacent index
njoin <- length(adj) # number of spatial joins
num <- sapply(nb, length) # number of neighbors of each grid


####################################################################################
####################################################################################
####################################################################################
# Model code in NIMBLE
distanceModelCode <- nimbleCode({
  ##### Priors #####
  for(i in 1:n_covar){
    beta[i] ~ dnorm(0, 0.001) # landscape var regression coefficient
  }
  muc ~ dunif(1, 3) # average cluster size
  sigma0 ~ dunif(0.1, 5) # sigma0
  p ~ dunif(0.1, 3) # model parameter for sigma ~ group_size
  
  # CAR priors
  tau ~ dgamma(0.001, 0.001) # precision param for spatial CAR
  weights[1:njoin] <- 1
  bspat[1:L] ~ dcar_normal(adj = adj[1:njoin], weights = weights[1:njoin], 
                           num = num[1:L], tau = tau, zero_mean = 0) # random spatial effect (intercept) intercept for each grid
  
  ##### Model #####
  
  ### Landscape abundance for each grid
  for(l in 1:L) {
    # Non-spatial abundance model
    log_lam[l] <- inprod(covar[l, 1:n_covar], beta[1:n_covar]) + bspat[l] # regression for log abundance (without intercept)
    lam[l] <- exp(log_lam[l]) # lambda for each grid
  }
  
  ### Abundance at each transect
  for (i in 1:I) {
    tr_lam[i] <- inprod(propM[i,1:L], lam[1:L])
  }
  
  ### Model prob group size class
  for(k in 1:K) {
    gs[k] <- (exp(-muc) * pow(muc, k)) / (factorial(k) * (1 - exp(-muc)))
  }
  
  ### Calculate sigma for each size class
  # sigma is function of group size class (k) bigger group size -> larger sigma[k] -> slower decay half normal detection function
  for(k in 1:K){
    log(sigma[k]) <- sigma0 + p * (k - 1)
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
    pi[k, 1] <- mn_cell[k, 1] * gs[k]
    for(j in 2:J){ # loop for distance class j >= 2
      # do the same but remove the cum prob of previous dist class j-1 instead
      mn_cell[k, j] <- (sqrt(2 * 3.1416) * sigma[k] / distBreaks[J]) * 
        (phi(distBreaks[j] / sigma[k]) - phi(distBreaks[j - 1] / sigma[k]))
      # p[k,j]
      pi[k, j] <- mn_cell[k, j] * gs[k]
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

constants <- list(I = ntran, J = dist_class_n, K = group_class_n, 
                  L = grid_n, n_covar = ncol(covar),
                  distBreaks = distBreaks,
                  adj = adj, num = num, njoin = njoin)
data <- list(y = y_matrix, covar = covar, propM = propM)
inits <- list(muc = runif(1, 1, 3), sigma0 = runif(1, 1, 2), p = runif(1, 0.5, 1), tau = runif(1,0.5,1),
              beta = rnorm(ncol(covar), 0,1), bspat = runif(grid_n, 0, 0.1))

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, 
                             constants = constants, 
                             data = data, 
                             inits = inits)

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = c("sigma", "p", "gs", "sigma0", "pi", 'lam', 'beta', 'tau', 'bspat'))
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

# Run the MCMC
#samples <- runMCMC(Cmcmc, niter = 60000, nburnin = 50000, thin = 5, nchains = 1)
###summary(samples)
#adjust to run MCMC niter = 6000, nburnin = 5000 #Apinya
samples <- runMCMC(Cmcmc, niter = 6, nburnin = 5, thin = 5, nchains = 1)
summary(samples)
################################ Post training ############################
# fitted detection function
sigma_1 <- median(samples[,'sigma[1]'])
sigma_2 <- median(samples[,'sigma[2]'])
sigma_3 <- median(samples[,'sigma[3]'])
sigma_4 <- median(samples[,'sigma[4]'])
p <- median(samples[,'p'])
p
sigma_5 <- median(samples[,'sigma[5]'])
sigma_6 <- median(samples[,'sigma[6]'])
# multinomial pi
pi_val <- apply(samples[,grep("^pi", colnames(samples))], 2, median)
pi <- matrix(pi_val, nrow=group_class_n, byrow = F)
pi
sum(pi)

### Visualize pi
plot(pi[1,], type='l', xlab='Distance Class', ylab='Group Size', main='Multinomial Prob')
lines(pi[2,], lty=2)
lines(pi[3,], lty=3)
legend('topright', legend=c('1','2','3'), title = 'Group size', lty=c(1,2,3))
#for BTG
lines(pi[4,], lty=4)
legend('topright', legend=c('1','2','3','4'), title = 'Group size', lty=c(1,2,3,4))

# calculate gx (only distance)
x_seq <- seq(0, dist_limit)
gx_1 <- exp((-1*x_seq^2)/(2*sigma_1^2))
gx_2 <- exp((-1*x_seq^2)/(2*sigma_2^2))
gx_3 <- exp((-1*x_seq^2)/(2*sigma_3^2))
hp <- hist(data_sub$P.dist, main='Histogram & Detection function', xlab='Distance')
####for BTG
gx_4 <- exp((-1*x_seq^2)/(2*sigma_4^2)) 
hp <- hist(data_sub$P.dist, main='Histogram & Detection function', xlab='Distance')
lines(x_seq, gx_1*max(hp$counts), col='red', lwd=2)
lines(x_seq, gx_2*max(hp$counts), col='green', lwd=2)
lines(x_seq, gx_3*max(hp$counts), col='blue', lwd=2)
####for BTG
lines(x_seq, gx_4*max(hp$counts), col='pink', lwd=4)
legend('topright', legend=c('1','2','3','4'), title = 'Group size', col=c('red', 'green', 'blue', 'pink'), lty=1)

3# Site-specific Lambda (Abundance)
lam <- samples[,grep("^lam", colnames(samples))]
tr_names <- dimnames(y_matrix)[[3]]
colnames(lam) <- 1:grid_n
boxplot(lam[,1:50], ylim=c(0,100))
boxplot(lam[,101:150], ylim=c(0,100))
abun_med <- apply(lam, 2, median, na.rm = TRUE)
hist(abun_med)

# Site-specific random effect
bspat <- samples[,grep("^bspat", colnames(samples))]
bspat_med <- apply(bspat, 2, median, na.rm = TRUE)

# Beta
beta <- samples[,grep("^beta", colnames(samples))]
colnames(beta) <- colnames(covar)
boxplot(beta)
abline(h=0, lty=3)

# Inject back to shp file
library(terra)
abun_med[abun_med>20]
poly$abun_med <- abun_med
poly$bspat_med <- bspat_med
plot(poly["abun_med"])
plot(poly["bspat_med"])
