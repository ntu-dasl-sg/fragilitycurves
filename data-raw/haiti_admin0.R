# This R script saves shapefiles of level 2 administrative units in Haiti:

library(rgdal)

# Read in shapefiles:
haiti_admin0 <- rgdal::readOGR("D:/Documents/Proj_Damage_Spatial_Corr/Data/hti_admbnda_adm0_cnigs_20181129.shp")

save(haiti_admin0, file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/haiti_admin0.RData")
