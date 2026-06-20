library(terra)
library(tidyverse)
library(raster)
NDVI_Nov2020<- raster("NDVI1km_NOV_2020.tif")
NDVI_Dec2020<- raster("NDVI1km_DEC_2020.tif")
NDVI_Jan2021<- raster("NDVI1km_JAN_2021.tif")
NDVI_Feb2021<- raster("NDVI1km_FEB_2021.tif")
NDVI_Mar2021<- raster("NDVI1km_MAR_2021.tif")

NDVI_Nov2020[NDVI_Nov2020 < 0] <- NA
NDVI_Dec2020[NDVI_Dec2020 < 0] <- NA
NDVI_Jan2021[NDVI_Jan2021 < 0] <- NA
NDVI_Feb2021[NDVI_Feb2021 < 0] <- NA
NDVI_Mar2021[NDVI_Mar2021 < 0] <- NA


rasterstack <- stack(NDVI_Nov2020,
                     NDVI_Dec2020,
                     NDVI_Jan2021,
                     NDVI_Feb2021,
                     NDVI_Mar2021)
plot(rasterstack)

NDVI_Mean <- calc(rasterstack, fun = mean, na.rm = TRUE)
NDVI_Std  <- calc(rasterstack, fun = sd, na.rm = TRUE)
NDVI_CV   <- NDVI_Std / NDVI_Mean

plot(NDVI_CV)

writeRaster(NDVI_CV, filename = "NDVI_CV_20250415.tif",overwrite=TRUE)

HKKSlope100m <- raster("D:/2024SpatialDistance/Covariate/HKKSlope100m.tif")

slope_Mean <-  terra::aggregate(HKKSlope100m, fact = 10, fun = "mean")
slope_Sd <-  terra::aggregate(HKKSlope100m, fact = 10, fun = "sd")
slope_var <- slope_Sd^2


slope_CV <- slope_Sd/slope_Mean
plot(slope_CV)

writeRaster(slope_var, filename = "slope_var.tif",overwrite=TRUE)
writeRaster(slope_CV, filename = "slope_CV.tif")


