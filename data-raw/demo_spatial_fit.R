# This R script fits a spatial ordinal model to the damage simulation data and saves the output:
library(MASS)
library(fragilitycurves) # With damage_simulation data and spatial_ordinal function included.

data(damage_simulation)

data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
data.subset.2 <- damage_simulation[damage_simulation$building_cat == 2, ]

frag.model.1 <- polr(CDF ~ logPGA, data = data.subset.1,
                     method = "probit", Hess = TRUE)
frag.model.2 <- polr(CDF ~ logPGA, data = data.subset.2,
                     method = "probit", Hess = TRUE)

# Set some reasonable upper and lower limits for parameters:

lower_lim <- rep(-Inf, 23); upper_lim <- rep(Inf, 23);

shared.range <- seq(0, 10, by = 1)
shared.vgm <- as.numeric(exp(-1.2313959))*(1-geoR::matern(shared.range, phi= as.numeric(exp(0.05)), kappa = 1)) + as.numeric(exp(-9.1425902))
shared.df <- data.frame("x" = shared.range, "y" = shared.vgm)

ggplot(data = shared.df, aes(x = x, y = y)) + geom_line() + labs(x = "km", y = expression(gamma(h))) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(a) Shared spatial field")

log_phi_max <- log(3); log_phi_min <- log(0.05);

# logphi: max of 3, min of 0.05.

log_slope1_max <- log(1.5*frag.model.1$coefficients);
log_slope1_min <- log(0.5*frag.model.2$coefficients)
log_slope2_max <- log(1.5*frag.model.1$coefficients);
log_slope2_min <- log(0.5*frag.model.2$coefficients)

cutoff.1.start <- frag.model.1$zeta
cutoff.2.start <- frag.model.2$zeta

# Reparameterising cut-offs to ensure increasing order in optimisation:
first_cutoff1 <-  cutoff.1.start[1]
first_cutoff2 <-  cutoff.2.start[1]

cutoff_factors <- function(cutoffs){
  temp <- rep(NA, length(cutoffs)-1)
  for (i in 2:length(cutoffs)){
    temp[i-1] <- cutoffs[i]-cutoffs[i-1]
  }
  return(temp)
}

cutoff_factors1 <- cutoff_factors(cutoff.1.start)
cutoff_factors2 <- cutoff_factors(cutoff.2.start)

cutoff11_max <- Inf; cutoff11_min <- -Inf
cutoff21_max <- Inf; cutoff21_min <- -Inf

factor_max <- 5*max(c(cutoff_factors1, cutoff_factors2));
factor_min <- 0.5*min(c(cutoff_factors1, cutoff_factors2))

lower_lim[1] <- log_phi_min;

lower_lim[10] <- log_slope1_min; lower_lim[11] <- log_slope2_min;

lower_lim[12:17] <- c(cutoff11_min, rep(factor_min, length(cutoff_factors1)));
lower_lim[18:23] <- c(cutoff21_min, rep(factor_min, length(cutoff_factors2)));

upper_lim[c(1, 4, 7)] <- log_phi_max;

upper_lim[10] <- log_slope1_max;  upper_lim[11] <- log_slope2_max;
upper_lim[12:17] <- c(cutoff11_max, rep(factor_max, length(cutoff_factors1)));
upper_lim[18:23] <- c(cutoff21_max, rep(factor_max, length(cutoff_factors2)));
upper_lim[3] <- -2;


# The following takes 8 minutes in a PC with characteristics: Intel(R) Xeon (R) W-2112 CPU Processor @ 3.60GHz; 32GB of RAM; Windows 10 64-bit:
temp.time <- proc.time()[3]
spatial_fit <- spatial_ordinal(frag.model.1, frag.model.2,
                               data.subset.1, data.subset.2,
                               lower.lim = lower_lim,
                               upper.lim = upper_lim)
time.taken <- proc.time()[3] - temp.time

demo_spatial_fit <- spatial_fit

save(demo_spatial_fit, file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/demo_spatial_fit.RData")
