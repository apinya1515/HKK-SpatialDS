library(raster)
library(tidyverse)

point<- shapefile("D:/WCS/GIS/HKK_Transect/point/32line_center.shp")

DEM<- raster("D:/WCS/GIS/HKK_Transect/resample/dem_resample.tif")

plot(DEM)
DEM_std <-  scale(DEM)
raster::extract(DEM_std,point)

DDF<- raster("D:/WCS/GIS/HKK_Transect/resample/ddf_resample.tif") %>% resample(DEM)
plot(DDF)
DDF_std <-  scale(DDF)
raster::extract(DDF_std,point)


DEF<- raster("D:/WCS/GIS/HKK_Transect/resample/def_resample.tif")%>% resample(DEM)
plot(DEF)
DEF_std <-  scale(DEF)
raster::extract(DEF_std,point)

MDF<- raster("D:/WCS/GIS/HKK_Transect/resample/mdf_resample.tif")%>% resample(DEM)
plot(MDF)
MDF_std <-  scale(MDF)
raster::extract(MDF_std,point)#extract 

STR<- raster("D:/WCS/GIS/HKK_Transect/Mainstream/MainStream_raster_dem1000_sum_nona_clip2.tif")%>% resample(DEM)
  
plot(STR)
STR_std <-  scale(STR)
str_extract <- raster::extract(STR_std,point)#extract 

stack_std<- stack(DEM_std,DDF_std, DEF_std, MDF_std,STR_std)

names(stack_std) <- c("ELE", "DDF", "DEF", "MDF","STR")


dataframe <- as.data.frame(point) %>% select(Line)

extract_dataframe <- bind_cols(raster::extract(stack_std, point))
  
extract_dataframe



# Distance sampling analysis in unmarked

#3 Importing, formatting, and summarizing data
library(unmarked)
SDdists <- read.csv("D:/DistanceUnmarked/SDdistdata.csv")
head(SDdists, 10)


yDat <- formatDistData(SDdists, distCol="distance",
                       transectNameCol="transect", dist.breaks=c(0, 10, 20, 30, 40,50,60,70,80,90,100,110,120))
yDat


covs <- extract_dataframe


umf <- unmarkedFrameDS(y=as.matrix(yDat), siteCovs=covs, survey="line",
                       dist.breaks=c(0, 10, 20, 30, 40,50,60,70,80,90,100,110,120), tlength=rep(32000, 32),
                       unitsIn="m")

summary(umf)
hist(umf, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)

#4 Model fitting
hn_Null <- distsamp(~1~1, umf)
hn_Null

hn_Null <- distsamp(~1~1, umf, output="density",
                    unitsOut="kmsq")
hn_Null
summary(hn_Null)
haz_Null <- distsamp(~1~1, umf, keyfun="hazard", output="density",
                     unitsOut="kmsq")
haz_Null
summary(haz_Null)

uni_Null <- distsamp(~1~1, umf, keyfun="uniform", output="density",
                     unitsOut="kmsq")
uni_Null
summary(uni_Null )

hn_Hab.Ht <- distsamp(~1~ELE+DDF+DEF+MDF, data= umf)
hn_Hab.Ht 

haz_Hab.Ht <- distsamp(~1~ELE+DDF+DEF+MDF, data= umf, keyfun= "hazard" ,output="density",unitsOut="kmsq")
uni_Hab.Ht <- distsamp(~1~ELE+DDF+DEF+MDF, data= umf, keyfun= "uniform" ,output="density",unitsOut="kmsq")

summary(hn_Hab.Ht)
summary(haz_Hab.Ht)
summary(uni_Hab.Ht)
hn_Hab.Ht <- distsamp(~1~ELE+DDF+DEF+MDF, data= umf, keyfun= "halfnorm")

#5 Manipulating results/ Back-transforming
names(haz_Null)
backTransform(haz_Null, type="state")
backTransform(haz_Null, type="det")
