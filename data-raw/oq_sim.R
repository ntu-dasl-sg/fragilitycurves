# This R script is for the selection of a portfolio of 150 buildings of each building type within a 2km x 2km grid in Port-au-Prince and the creation of the dataframe containing its OpenQuake logPGA simulations:

library(sp)
library(raster)

load(file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/logPGA_df_OQ.RData")
# logPGA_df.

head(logPGA_df[, 1:10])
events_affected <- sum(!(colnames(logPGA_df) %in% c("lon", "lat", "building_cat")))
# 156491.

pap.pts <- SpatialPoints(logPGA_df[, c("lon", "lat")], proj4string = crs(pap_shp))

# 2km x 2km:
study_area <- Polygon(rbind(c(-72.34, 18.535), c(-72.34, 18.555),
                            c(-72.32, 18.555), c(-72.32, 18.535)))
study_area <- Polygons(list(study_area),1)
study_area <- SpatialPolygons(list(study_area))
crs(study_area) <- crs(pap.pts)

# Keep points in area only:
study_ind <- over(pap.pts, study_area)

study_data <- logPGA_df[!is.na(study_ind), ]

set.seed(1)
oq_sim <- study_data  %>% group_by(building_cat) %>% sample_n(150)

save(oq_sim, file = "F:/Documents/RPackages/fragilitycurves/data/oq_sim.RData")
