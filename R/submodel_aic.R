#' Compute Akaike Information Criterion (AIC) values for the Damage-spatial model and its Non-spatial and IM-spatial submodels.
#'
#' @param data.1 The damage dataset of Building Category 1 used in the fit.
#' @param data.2 The damage dataset of Building Category 2 used in the fit.
#' @param model.fit The output of a \code{spatial_ordinal} model fit.
#' @return A vector containing the AIC values for the Non-spatial and IM-spatial submodels as well as the Damage-spatial model.
#' @importFrom stats pnorm
#' @export
#'
#' @examples
#' data(demo_spatial_fit)
#' data(damage_simulation)
#'
#' data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
#' data.subset.2 <- damage_simulation[damage_simulation$building_cat == 2, ]
#'
#' # The following takes 10 seconds on a PC with characteristics:
#' # Intel(R) Xeon (R) W-2112 CPU Processor @ 3.60GHz; 32GB of RAM;
#' aic_val <- submodel_aic(data.1 = data.subset.1,
#'                          data.2 = data.subset.2, model.fit = demo_spatial_fit)

submodel_aic <- function(data.1, data.2, model.fit){

  data <- rbind(data.1, data.2)
  CDF_breaks <- sort(unique(data.1$CDF), decreasing = FALSE)

  dist.mat <-as.matrix(stats::dist(data[, c('Easting', 'Northing')], method = "euclidean", diag = TRUE))
  dist.mat <- dist.mat/1000 # Work in km instead of m.
  dist.mat.1 <-as.matrix(stats::dist(data.1[, c('Easting', 'Northing')], method = "euclidean", diag = TRUE))
  dist.mat.2 <-as.matrix(stats::dist(data.2[, c('Easting', 'Northing')], method = "euclidean", diag = TRUE))
  dist.mat.1 <- dist.mat.1/1000 # Work in km instead of m.
  dist.mat.2 <- dist.mat.2/1000 # Work in km instead of m.

  no_states <- nlevels(data.1$CDF)
  damage_ind_1 <- matrix(NA, nrow = nrow(data.1), ncol = no_states)
  damage_ind_2 <- matrix(NA, nrow = nrow(data.2), ncol = no_states)
  for (i in 1:no_states){
    damage_ind_1[, i] <- as.numeric(data.1$CDF == CDF_breaks[i])
    damage_ind_2[, i] <- as.numeric(data.2$CDF == CDF_breaks[i])
  }
  # Check:
  # which(rowSums(damage_ind_1)!=1)
  # which(rowSums(damage_ind_2)!=1)

  damage_lower_1 <- damage_ind_1[, 2:no_states]; damage_lower_2 <- damage_ind_2[, 2:no_states]
  damage_ind_1 <- damage_ind_1[, 1:(no_states-1)]; damage_ind_2 <- damage_ind_2[, 1:(no_states-1)]

  damage_state_1 <- as.numeric(data.1$CDF) # Previously ordered.
  damage_state_2 <- as.numeric(data.2$CDF) # Previously ordered.

  input_data <- list(
    D = dist.mat,
    x1 = data.1$logPGA,
    x2 = data.2$logPGA,
    damage_ind_1 = damage_ind_1,
    damage_lower_1 = damage_lower_1,
    damage_state_1 = damage_state_1,
    damage_ind_2 = damage_ind_2,
    damage_lower_2 = damage_lower_2,
    damage_state_2 = damage_state_2,
    no_states = no_states,
    log_kappa = 0 # Fix smoothness parameter to INLA default to separate noise and smoothness.
  )

  pgaspat_par <- list(
    c_factor1 = model.fit$par[names(model.fit$par) == "c_factor1"],
    c_factor2 = model.fit$par[names(model.fit$par) == "c_factor2"],
    log_phi = model.fit$par["log_phi"],
    log_sigma_2 =  model.fit$par["log_sigma_2"],
    log_tau_2 =  model.fit$par["log_tau_2"],
    log_tau1_2 =  log(exp(model.fit$par["log_tau1_2"]) + exp(model.fit$par["log_sigma1_2"]) + exp(model.fit$par["log_tau_2"])),
    log_tau2_2 =  log(exp(model.fit$par["log_tau2_2"]) + exp(model.fit$par["log_sigma2_2"]) + exp(model.fit$par["log_tau_2"])),
    field = model.fit$field,
    log_slope1 =  model.fit$par["log_slope1"],
    log_slope2 =  model.fit$par["log_slope2"])

  pgaspat_obj <- MakeADFun(
    data = input_data,
    parameters = pgaspat_par,
    random = 'field',
    DLL = "pgaspatial_newmarginals_shared")

  # Test functions:
  pga_spat_opt <- pgaspat_obj$fn(pgaspat_obj$par)

  nonspat_par <- list(
    c_factor1 = model.fit$par[names(model.fit$par) == "c_factor1"],
    c_factor2 = model.fit$par[names(model.fit$par) == "c_factor2"],
    log_tau_2 =  model.fit$par["log_tau_2"],
    log_tau1_2 =  log(exp(model.fit$par["log_tau1_2"]) + exp(model.fit$par["log_sigma1_2"]) + exp(model.fit$par["log_tau_2"])),
    log_tau2_2 =  log(exp(model.fit$par["log_tau2_2"]) + exp(model.fit$par["log_sigma2_2"]) + exp(model.fit$par["log_tau_2"])),
    field = model.fit$field,
    log_slope1 =  model.fit$par["log_slope1"],
    log_slope2 =  model.fit$par["log_slope2"])

  nonspat_obj <- MakeADFun(
    data = input_data,
    parameters = nonspat_par,
    random = 'field',
    DLL = "nonspatial_newmarginals_shared")

  # Test functions:
  non_spat_opt <-nonspat_obj$fn(nonspat_obj$par)

  nonspat_aic <- 2*non_spat_opt + 2*(length(model.fit$par)-(2*2+2))

  pga_spat_aic <- 2*pga_spat_opt + 2*(length(model.fit$par)-2*2)

  damage_spat_aic <- 2*model.fit$objective + 2*length(model.fit$par)

  temp_vec <- c(nonspat_aic, pga_spat_aic, damage_spat_aic)
  names(temp_vec) <- c("Non-spatial", "IM-spatial", "Damage-spatial")

  return(temp_vec)

}
