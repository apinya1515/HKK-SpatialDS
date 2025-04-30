#### This code consider group size as categorical

library(nimble)
library(dplyr)

# Transect data (table 1)
data_tr <- read.table('D:\\UserData\\Dropbox\\@KU_pop_course\\day1_exercise\\line_data.txt', 
                       sep='\t', header=T)

##############################################################################
### Input parameters
species <- 'SBR'
dist_limit <- 90 # max observed distance
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
data_sub <- data_tr %>% filter(Species==species) # Sambar deer

group_max <- max(data_sub$Gz.sz)   # max group size
group_class_size <- seq(1,group_max)  # group size in each class (1 increment)
group_class_n <- length(group_class_size) # number of group size classes

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
  print(paste0('data ',i, ' distance:',  round(samp$P.dist,2),' in grsz:', k, ' dist:', j, ' transect:', tr))
  y_matrix[k, j, tr] = y_matrix[k, j, tr] + 1 # add observation into group class - dist class - transect
}
print(y_matrix)

sum(y_matrix) # check the number of observations
nrow(data_sub) # check the number of observations

####################################################################################
# Model code in NIMBLE
distanceModelCode <- nimbleCode({
  ##### Priors #####
  muc ~ dunif(0, 30) # average cluster size (lamdas)
  sigma0 ~ dunif(0, 10) # sigma0
  p ~ dunif(0, 10) # model parameter for sigma ~ group_size
  ### Priors for expected number of clusters for each transect
  for(i in 1:I){
    lam[i] ~ dunif(0, 100) # lambda for abundance of clusters for each transect
  }
  
  ##### Model #####
  
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
    # pi[k,j] = probability that the obs is in class k and j
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
        mu[k, j, i] <- lam[i] * pi[k, j] # probability of class k and j * abundance of transect i
        y[k, j, i] ~ dpois(mu[k, j, i])
      }
    }
  }
})

####################################################################################

# Initial values
inits <- list(muc = runif(1, 1, 2),
              sigma0 = runif(1, 1, 5),
              p = runif(1, 0.01, 1))

# Constants
constants <- list(I = ntran,  # Number of sites (transect)
                  J = dist_class_n,  # Number of distance classes
                  K = group_class_n, # Number of group size classes
                  X = dist_class_mid,
                  distBreaks = distBreaks)
# Data
data <- list(y = y_matrix)

# Build the model
distanceModel <- nimbleModel(code = distanceModelCode, 
                             constants = constants, 
                             data = data, 
                             inits = inits)

# Configure and compile the MCMC
mcmcConf <- configureMCMC(distanceModel, monitors = c("sigma", "muc", "p", "gs", "sigma0", "pi", 'lam'))
distanceMCMC <- buildMCMC(mcmcConf)
CdistanceModel <- compileNimble(distanceModel)
Cmcmc <- compileNimble(distanceMCMC, project = distanceModel)

# Run the MCMC
samples <- runMCMC(Cmcmc, niter = 50000, nburnin = 40000, thin = 5, nchains = 1)
summary(samples)

# fitted detection function
sigma_1 <- median(samples[,'sigma[1]'])
sigma_2 <- median(samples[,'sigma[2]'])
sigma_3 <- median(samples[,'sigma[3]'])
p <- median(samples[,'p'])

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

# calculate gx (only distance)
x_seq <- seq(0, dist_limit)
gx_1 <- exp((-1*x_seq^2)/(2*sigma_1^2))
gx_2 <- exp((-1*x_seq^2)/(2*sigma_2^2))
gx_3 <- exp((-1*x_seq^2)/(2*sigma_3^2))
hp <- hist(data_sub$P.dist, main='Histogram & Detection function', xlab='Distance')
lines(x_seq, gx_1*max(hp$counts), col='red', lwd=2)
lines(x_seq, gx_2*max(hp$counts), col='green', lwd=2)
lines(x_seq, gx_3*max(hp$counts), col='blue', lwd=2)
legend('topright', legend=c('1','2','3'), title = 'Group size', col=c('red', 'green', 'blue'), lty=1)

# Site-specific Lambda (Abundance)
lam <- samples[,grep("^lam", colnames(samples))]
tr_names <- dimnames(y_matrix)[[3]]
colnames(lam) <- tr_names
boxplot(lam)

