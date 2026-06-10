##############################################################################

### Get data from files
tran_id <- sort(unique(data_tr$Tr.no)) # name of transect
tran_n <- length(tran_id) # number of transect
grid_id <- data_land$grid_id # grid id
grid_n <- nrow(data_land) # number of grid
covar <- as.data.frame(data_land[,-1] %>% scale()) # select only covariate columns and rescale

### Create observation matrix (matrix=transect, row=cluster size class ,column=distance class)
dist_width <- dist_limit/dist_class_n # interval between distance class
distBreaks = seq(dist_width, dist_class_n*dist_width, by=dist_width) # upper breaks of distance classes

# subset data by species
data_sub <- data_tr %>% filter(Species==species)

gs_max <- max(data_sub$Gz.sz)   # max group size
gs_class_n <- length(gsBreaks) # number og gs classes
# Ensure the last break in gsBreaks goes up to the maximum observed group size (gs_max)
# to avoid truncation error and cover all observed group sizes in the final category.
gsBreaks[gs_class_n] <- gs_max

#### Create zeros-array for the data
y_matrix <- array(0, dim=c(gs_class_n, dist_class_n, tran_n), 
                  dimnames=list(1:gs_class_n, distBreaks, tran_id)) # 3 dimensions = group size, distance, transect

### Populate y data for target species ... by +1 to y_matrix for each data point that fall into the cell
# y_matrix[group_size, distance, transect]
for(i in 1:nrow(data_sub)){
  samp <- data_sub[i,] # get data
  tr <- which(tran_id ==  samp$Tr.no) # transect index
  j <- sum(distBreaks < samp$P.dist) + 1 # distance class index
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
lenM <- matrix(0, nrow=tran_n, ncol=grid_n) # create zero matrix row=#of trasect col=#of grid
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
# Create a neighbors list (Queen's contiguity = 8 neighbors)
nb <- poly2nb(poly, queen = TRUE, row.names = poly$ID)
adj <- unlist(nb) # adjacency index (grid id of each neighbor)
njoin <- length(adj) # number of spatial joins
num <- sapply(nb, length) # number of neighbors of each grid 
### (sum of num must match with the length of adj)

## pre-calculate factorial
# lgamma(x) is equivalent to the natural-log of the factorial = log((x-1)!)
# So, lgamma(x+1) = log(x!)
# Calculate up to gs_max to prevent index-out-of-bounds in the model loop
logFactorial = lgamma((1:gs_max) + 1)

# Define a water mask: 0 if grid cell is mostly water (WA > 0.5), 1 otherwise
water_mask <- ifelse(data_land_orig$WA > 0.5, 0, 1)

