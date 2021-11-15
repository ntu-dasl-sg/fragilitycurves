# This R script simulates sample datasets for spatially correlated damage for Haiti:

rm(list = ls())
library(gstat)
library(raster)
# library(malariaAtlas) # Loads ggplot2.
library(dplyr)
library(ggplot2)
library(lemon)

# Read in fitted spatial model parameters:
load("D:/Documents/Proj_Damage_Spatial_Corr/Data/haiti_feb_2_mega.RData")
rm(list = ls()[!(ls()%in%c("par.est", "data.path"))])
par.est

# load(file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/logPGA_df_OQ.RData") # logPGA_df

# Use simulated mean PGA for 2010 event:
mean.PGA.raster <- raster("D:/Documents/Proj_Damage_Spatial_Corr/Data/Haiti_shaking_Raster_v2.grd")
mean.PGA.raster <- mean.PGA.raster*100 # Change to percentage of gravity.
crs(mean.PGA.raster) <- CRS("+proj=longlat")
mean.PGA.raster <- projectRaster(from = mean.PGA.raster,
                                 crs = CRS("+proj=utm +zone=18 ellps=WGS84"))

# Crop raster by Haiti outline:
haiti_admin0 <- rgdal::readOGR("D:/Documents/Proj_Damage_Spatial_Corr/Data/hti_admbnda_adm0_cnigs_20181129.shp")
study_shp <- sp::spTransform(haiti_admin0, CRS("+proj=utm +zone=18 ellps=WGS84"))
cropped.raster <- mask(mean.PGA.raster, study_shp)

# Use Mega Subset 1's estimates:
mega_subset <- 1

field.phi <- exp(par.est$log_phi[mega_subset]); field.sigma2 <- exp(par.est$log_sigma_2[mega_subset]); field.tau2 <- exp(par.est$log_tau_2[mega_subset]);
field1.phi <- exp(par.est$log_phi1[mega_subset]); field1.sigma2 <- exp(par.est$log_sigma1_2[mega_subset]); field1.tau2 <- exp(par.est$log_tau1_2[mega_subset]);
field2.phi <- exp(par.est$log_phi2[mega_subset]); field2.sigma2 <- exp(par.est$log_sigma2_2[mega_subset]); field2.tau2 <- exp(par.est$log_tau2_2[mega_subset]);
beta1 <- exp(par.est$log_slope1[mega_subset]);
beta2 <- exp(par.est$log_slope2[mega_subset]);

cutoffs_factor1 <- par.est[mega_subset, names(par.est) == "c_factor1"]
cutoffs1 <- as.numeric(cutoffs_factor1)
names(cutoffs1) <- rep("cutoffs1", length(cutoffs1))

no_states <- length(cutoffs1) + 1

for (i in 2:(no_states - 1)){
  cutoffs1[i] <-  cutoffs1[i]*prod(cutoffs_factor1[1:(i-1)])
}

cutoffs_factor2 <- par.est[mega_subset, names(par.est) == "c_factor2"]
cutoffs2 <- as.numeric(cutoffs_factor2)
names(cutoffs2) <- rep("cutoffs2", length(cutoffs2))

for (i in 2:(no_states - 1)){
  cutoffs2[i] <-  cutoffs2[i]*prod(cutoffs_factor2[1:(i-1)])
}

# Simulate fields:

UTM.pts <- rasterToPoints(cropped.raster)
colnames(UTM.pts) <- c("Easting", "Northing", "Val")
UTM.pts <- as.data.frame(UTM.pts)

# Use a balanced sample of PGA values:
UTM.pts$PGA_bin <- cut(UTM.pts$Val, 6)

table(UTM.pts$PGA_bin)

UTM.pts <- rbind(UTM.pts, UTM.pts)

UTM.pts$building_cat <- rep(c(1, 2), each = nrow(UTM.pts)/2)

nrow(UTM.pts)

head(UTM.pts)

tail(UTM.pts)

set.seed(1)

UTM.pts <- UTM.pts %>% group_by(PGA_bin, building_cat) %>% sample_n(100, replace = TRUE) # Replace = TRUE to allow for nugget estimation.

UTM.pts[, c("Easting", "Northing")] <- UTM.pts[, c("Easting", "Northing")]/1000 # Parameters relate to km for distance not m.

nrow(UTM.pts)

# a. Shared field:

# define the gstat object (spatial model)
g.dummy <- gstat(formula=z~1, locations=~Easting+Northing, dummy=T, beta=0, model=vgm(psill=field.sigma2,range=field.phi,nugget=field.tau2,kappa=1,model="Mat"))

# Make simulations based on the stat object
set.seed(2)
temp.time <- proc.time()[3]
field.sim <- predict(g.dummy, newdata=UTM.pts, nsim=1)
time.taken.2 <- proc.time()[3] - temp.time # 1.84s.

field.loc.mean <- rowMeans(field.sim@data)
hist(field.loc.mean)

# b. Field for Building Cat 1:

# Define the gstat object (spatial model)
g.dummy1 <- gstat(formula=z~1, locations=~Easting+Northing, dummy=T, beta=0, model=vgm(psill=field1.sigma2,range=field1.phi,nugget=0,kappa=1,model="Mat"))

set.seed(3)
temp.time <- proc.time()[3]
field1.sim <- predict(g.dummy1, newdata=UTM.pts[UTM.pts$building_cat == 1, ], nsim=1)
time.taken.3 <- proc.time()[3] - temp.time
# About 0.31s.

# c. Field for Building Cat 2:

# Define the gstat object (spatial model)
g.dummy2 <- gstat(formula=z~1, locations=~Easting+Northing, dummy=T, beta=0, model=vgm(psill=field2.sigma2,range=field2.phi,nugget=0,kappa=1,model="Mat"))

set.seed(4)
temp.time <- proc.time()[3]
field2.sim <- predict(g.dummy2, newdata=UTM.pts[UTM.pts$building_cat == 2, ], nsim=1)
time.taken.3 <- proc.time()[3] - temp.time
# About 0.39s.

# 5. Compute latent variables and prob of collapse:

# Convert NA logPGA values to -Inf:
UTM.pts$logPGA <- log(UTM.pts$Val)
logPGA_sim <- UTM.pts$logPGA
logPGA_sim[is.na(logPGA_sim)] <- -Inf
logPGA_sim

# Generate random error terms to be added to the non-spatial and PGA-spatial latent variables only.
set.seed(5)

iid.error3 <- matrix(rnorm(n = nrow(UTM.pts), mean = 0, sd = sqrt(field.tau2 + field.sigma2)), nrow = nrow(UTM.pts), ncol = 1)

# Latent variables:

nonspat.lv1 <- beta1*logPGA_sim[UTM.pts$building_cat == 1]

nonspat.lv2 <- beta2*logPGA_sim[UTM.pts$building_cat == 2]

sharedspat.lv1 <- nonspat.lv1 + beta1*field.sim$sim1[UTM.pts$building_cat == 1]

sharedspat.lv2 <- nonspat.lv2 + beta2*field.sim$sim1[UTM.pts$building_cat == 2]

damagespat.lv1 <- sharedspat.lv1 + field1.sim$sim1

damagespat.lv2 <- sharedspat.lv2 + field2.sim$sim1

UTM.pts$lv <- NA
UTM.pts$lv[UTM.pts$building_cat == 1] <- damagespat.lv1
UTM.pts$lv[UTM.pts$building_cat == 2] <- damagespat.lv2

# Rescale Easting and Northing back:

UTM.pts[, c("Easting", "Northing")] <-UTM.pts[, c("Easting", "Northing")]*1000

ggplot() + geom_point(data = UTM.pts[UTM.pts$building_cat == 1, ], aes(x = Easting, y = Northing, color = lv)) + geom_polygon(data = study_shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

ggplot() + geom_point(data = UTM.pts[UTM.pts$building_cat == 2, ], aes(x = Easting, y = Northing, color = lv)) + geom_polygon(data = study_shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)


damage.sim <- function(cutoffs, lv, sd_val){
  # cutoffs <- cutoffs1; lv <- damagespat.lv1[, 1]; sd_val <- sqrt(field1.tau2)
  simulated_damage <- rep(NA, length(lv))
  temp <- matrix(0, nrow = length(lv), ncol = no_states) # Compute the damage probabilities per state.
  for (j in 1:(no_states-2)){ # Compute the probability of each damage state for all locations.
    temp[, j+1] <-   pnorm((cutoffs[j+1] - lv)/sd_val) - pnorm((cutoffs[j] - lv)/sd_val)
  }
  temp[, no_states] <- 1 - pnorm((cutoffs[no_states-1] - lv)/sd_val)
  temp[, 1] <- 1 - rowSums(temp)
  # Compute the mean replacement cost as a weighted average:

  for (i in 1:length(simulated_damage)){
    simulated_damage[i] <- sample(1:no_states, size = 1, prob = temp[i,])
  }

  return(simulated_damage)
}

set.seed(6)
damagespat.sim1 <- damage.sim(cutoffs1, damagespat.lv1, sqrt(field1.tau2))

set.seed(7)
damagespat.sim2 <- damage.sim(cutoffs2, damagespat.lv2, sqrt(field2.tau2))

UTM.pts$PGA <- exp(UTM.pts$logPGA)
UTM.pts$CDF <- c(damagespat.sim1, damagespat.sim2)

damage.data <- read.csv("D:/Documents/Proj_Damage_Spatial_Corr/Data/new_Haiti_Data.csv")

CDF_breaks <- sort(unique(damage.data$CDF), decreasing = FALSE)

UTM.pts$CDF <- ordered(UTM.pts$CDF, levels = 1:no_states, labels = as.character(CDF_breaks))

damage_simulation <- UTM.pts[, c("building_cat", "CDF", "logPGA", "PGA", "Easting", "Northing")]

# Create a balanced dataset of each building category:

damage_simulation_1 <- damage_simulation[damage_simulation$building_cat == "1", ]

table(damage_simulation_1$CDF)

set.seed(8)
damage_simulation_1 <- damage_simulation_1 %>% group_by(CDF) %>% sample_n(25)

damage_simulation_2 <- damage_simulation[damage_simulation$building_cat == "2", ]

table(damage_simulation_2$CDF)

set.seed(9)
damage_simulation_2 <- damage_simulation_2 %>% group_by(CDF) %>% sample_n(25)

damage_simulation <- rbind(damage_simulation_1, damage_simulation_2)

save(damage_simulation, file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/damage_simulation.RData")

# Check:

load(file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/damage_simulation.RData")

library(MASS)

data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
data.subset.2 <- damage_simulation[damage_simulation$building_cat == 2, ]

frag.model.1 <- polr(CDF ~ logPGA, data = data.subset.1,
                     method = "probit", Hess = TRUE)
frag.model.2 <- polr(CDF ~ logPGA, data = data.subset.2,
                     method = "probit", Hess = TRUE)

frag.model.1
frag.model.2

damage_simulation$building_cat <- as.character(damage_simulation$building_cat)

plot1 <- ggplot() + geom_point(data = damage_simulation_1, aes(x = Easting, y = Northing, color = CDF)) + geom_polygon(data = study_shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) + ggtitle("Building Category 1")

plot2 <- ggplot() + geom_point(data = damage_simulation_2, aes(x = Easting, y = Northing, color = CDF)) + geom_polygon(data = study_shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA) + ggtitle("Building Category 2")

grid_arrange_shared_legend(plot1, plot2, ncol = 2)
