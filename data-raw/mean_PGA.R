# This R script saves a simulated Haiti shaking map as .RData:

library(raster)
library(rgdal)
library(sp)

# Read in mean PGA surface:
mean.PGA.raster <- raster("D:/Documents/Proj_Damage_Spatial_Corr/Data/Haiti_shaking_Raster_v2.grd")
mean.PGA.raster <- mean.PGA.raster*100 # Change to percentage of gravity.
crs(mean.PGA.raster) <- CRS("+proj=longlat")
mean.PGA.raster <- projectRaster(from = mean.PGA.raster,
                                  crs = CRS("+proj=utm +zone=18 ellps=WGS84"))

mean_PGA <- mean.PGA.raster

save(mean_PGA, file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/mean_PGA.RData")
