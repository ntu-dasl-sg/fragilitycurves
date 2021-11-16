#' Fit a spatial ordinal model based on damage data for two building categories.
#'
#' @param model.1 An non-spatial ordinal model fitted for \code{data.1} (Building Category 1) using \code{polr} from the R package \code{MASS}.
#' @param model.2 An non-spatial ordinal model fitted for \code{data.2} (Building Category 2) using \code{polr} from the R package \code{MASS}.
#' @param data.1 A dataframe with the columns \code{CDF} (ordered factor), \code{logPGA} (numeric), \code{PGA} (numeric), \code{Easting} (numeric) and \code{Northing} (numeric) for Building Category 1.
#' @param data.2 A dataframe with the columns \code{CDF} (ordered factor), \code{logPGA} (numeric), \code{PGA} (numeric), \code{Easting} (numeric) and \code{Northing} (numeric) for Building Category 2.
#' @param start A list of starting values for the model parameters: \code{c_factor1} (first cut-off value and cut-off ratios for Building Category 1; this is a vector with one less component than the number of damage grades), \code{c_factor2} (first cut-off value and cut-off ratios for Building Category 2; this is a vector with one less component than the number of damage grades), \code{log_phi} (log-transformed range parameter for the shared spatial field), \code{log_sigma_2} (log-transformed partial sill parameter for the shared spatial field), \code{log_phi1} (log-transformed range parameter for the Building Category 1 spatial field), \code{log_tau1_2} (log-transformed nugget parameter for the Building Category 1 spatial field), \code{log_sigma1_2} (log-transformed partial sill parameter for the Building Category 1 spatial field), \code{log_phi2} (log-transformed range parameter for the Building Category 2 spatial field), \code{log_tau2_2} (log-transformed nugget parameter for the Building Category 2 spatial field), \code{log_sigma2_2} (log-transformed partial sill parameter for the Building Category 2 spatial field), \code{log_slope1} (log-transformed slope parameter for Building Category 1) and \code{log_slope2} (log-transformed slope parameter for Building Category 2).
#' @param lower.lim A vector of lower bounds for the parameters (in this order): \code{log_phi}, \code{log_sigma_2}, \code{log_tau_2}, \code{log_phi1}, \code{log_tau1_2}, \code{log_sigma1_2}, \code{log_phi2},  \code{log_tau2_2}, \code{log_sigma2_2}, \code{log_slope1}, \code{log_slope2}, \code{c_factor1}, \code{c_factor1}, \code{c_factor1}, \code{c_factor1}, \code{c_factor1}, \code{c_factor1}, \code{c_factor2}, \code{c_factor2}, \code{c_factor2}, \code{c_factor2}, \code{c_factor2}, \code{c_factor2}.
#' @param upper.lim A vector of upper bounds for the parameters with the parameter order corresponding to that in \code{lower.lim}.
#' @return A list of the model fit results including the parameter estimates, the number of iterations taken and the \code{nlminb} convergence code. This is the output of the \code{nlminb} function applied to minimise the negative log-likelihood of the spatial ordinal model together with the field estimates and the estimated parameter standard errors.
#' @importFrom TMB MakeADFun
#' @importFrom TMB sdreport
#' @importFrom stats nlminb
#' @export
#'
#' @examples
#' library(MASS)
#' data(damage_simulation)
#'
#' data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
#' data.subset.2 <- damage_simulation[damage_simulation$building_cat == 2, ]
#'
#' frag.model.1 <- polr(CDF ~ logPGA, data = data.subset.1,
#'                     method = "probit", Hess = TRUE)
#' frag.model.2 <- polr(CDF ~ logPGA, data = data.subset.2,
#'                     method = "probit", Hess = TRUE)
#'
#' # Set some reasonable upper and lower limits for parameters:
#' lower_lim <- rep(-Inf, 23); upper_lim <- rep(Inf, 23);
#'
#' log_phi_max <- log(3); log_phi_min <- log(0.05);
#'
#' log_slope1_max <- log(1.5*frag.model.1$coefficients);
#' log_slope1_min <- log(0.5*frag.model.2$coefficients)
#' log_slope2_max <- log(1.5*frag.model.1$coefficients);
#' log_slope2_min <- log(0.5*frag.model.2$coefficients)
#'
#' cutoff.1.start <- frag.model.1$zeta
#' cutoff.2.start <- frag.model.2$zeta
#'
#' # Reparameterising cut-offs to ensure increasing order in optimisation:
#' first_cutoff1 <-  cutoff.1.start[1]
#' first_cutoff2 <-  cutoff.2.start[1]
#'
#' cutoff_factors <- function(cutoffs){
#'  temp <- rep(NA, length(cutoffs)-1)
#'  for (i in 2:length(cutoffs)){
#'    temp[i-1] <- cutoffs[i]-cutoffs[i-1]
#'  }
#'  return(temp)
#' }
#'
#' cutoff_factors1 <- cutoff_factors(cutoff.1.start)
#' cutoff_factors2 <- cutoff_factors(cutoff.2.start)
#'
#' factor_max <- 5*max(c(cutoff_factors1, cutoff_factors2));
#' factor_min <- 0.5*min(c(cutoff_factors1, cutoff_factors2))
#'
#' lower_lim[1] <- log_phi_min;
#'
#' lower_lim[10] <- log_slope1_min; lower_lim[11] <- log_slope2_min;
#'
#' lower_lim[12:17] <- c(-Inf, rep(factor_min, length(cutoff_factors1)));
#' lower_lim[18:23] <- c(-Inf, rep(factor_min, length(cutoff_factors2)));
#'
#' upper_lim[10] <- log_slope1_max;  upper_lim[11] <- log_slope2_max;
#' upper_lim[12:17] <- c(Inf, rep(factor_max, length(cutoff_factors1)));
#' upper_lim[18:23] <- c(Inf, rep(factor_max, length(cutoff_factors2)));
#' upper_lim[3] <- -2;
#'
#' # The following takes 8 minutes on a PC with characteristics:
#' # Intel(R) Xeon (R) W-2112 CPU Processor @ 3.60GHz; 32GB of RAM; Windows 10 64-bit:
#' #temp.time <- proc.time()[3]
#' #spatial_fit <- spatial_ordinal(frag.model.1, frag.model.2,
#' #                              data.subset.1, data.subset.2,
#' #                              lower.lim = lower_lim,
#' #                             upper.lim = upper_lim)
#' #time.taken <- proc.time()[3] - temp.time

spatial_ordinal <- function(model.1, model.2, data.1, data.2, start = NULL, lower.lim = rep(-Inf, 23), upper.lim = rep(Inf, 23)){

  if(is.null(start)){

    slope.1.start <- max(model.1$coefficients, 0.01)
    slope.2.start <- max(model.2$coefficients, 0.01)
    cutoff.1.start <- model.1$zeta
    cutoff.2.start <- model.2$zeta

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

    parameters <- list(
      c_factor1 = c(first_cutoff1, cutoff_factors1),
      c_factor2 = c(first_cutoff2, cutoff_factors2),
      log_phi = 1.6532815,
      log_sigma_2 = -1.2313959,
      log_tau_2 = -9.1425902,
      log_phi1 = -1.1385036,
      log_tau1_2 = -0.7311687,
      log_sigma1_2 = -2.8401772,
      log_phi2 = -0.7887266,
      log_tau2_2 = -0.8261514,
      log_sigma2_2 = -1.2392615,
      field = rep(0, nrow(data.subset.1) + nrow(data.subset.2)),
      field1 = rep(0, nrow(data.subset.1)),
      field2 = rep(0, nrow(data.subset.2)),
      log_slope1 = log(slope.1.start),
      log_slope2 = log(slope.2.start)
    )

  }else{
    parameters <- start; parameters$log_tau_2 <- -9.1425902;
    parameters$field <- rep(0, nrow(data.subset.1) + nrow(data.subset.2));
    parameters$field <- rep(0, nrow(data.subset.1) + nrow(data.subset.2));
    parameters$field <- rep(0, nrow(data.subset.1) + nrow(data.subset.2));
    parameters$log_slope1 <- log(slope.1.start);
    parameters$log_slope2 <- log(slope.2.start);
  }

  data.subset <- rbind(data.1, data.2)

  dist.mat <-as.matrix(stats::dist(data.subset[, c('Easting', 'Northing')], method = "euclidean", diag = TRUE))
  dist.mat <- dist.mat/1000 # Work in km instead of m.
  dist.mat.1 <-as.matrix(stats::dist(data.1[, c('Easting', 'Northing')], method = "euclidean", diag = TRUE))
  dist.mat.2 <-as.matrix(stats::dist(data.2[, c('Easting', 'Northing')], method = "euclidean", diag = TRUE))
  dist.mat.1 <- dist.mat.1/1000 # Work in km instead of m.
  dist.mat.2 <- dist.mat.2/1000 # Work in km instead of m.

  no_states <- nlevels(data.1$CDF)
  CDF_breaks <- as.character(levels(data.1$CDF))

  damage_ind_1 <- matrix(NA, nrow = nrow(data.1), ncol = no_states)
  damage_ind_2 <- matrix(NA, nrow = nrow(data.2), ncol = no_states)
  for (i in 1:no_states){
    damage_ind_1[, i] <- as.numeric(data.1$CDF == CDF_breaks[i])
    damage_ind_2[, i] <- as.numeric(data.2$CDF == CDF_breaks[i])
  }

  damage_lower_1 <- damage_ind_1[, 2:no_states]; damage_lower_2 <- damage_ind_2[, 2:no_states]
  damage_ind_1 <- damage_ind_1[, 1:(no_states-1)]; damage_ind_2 <- damage_ind_2[, 1:(no_states-1)]

  damage_state_1 <- as.numeric(data.1$CDF) # Previously ordered.
  damage_state_2 <- as.numeric(data.2$CDF) # Previously ordered.

  input_data <- list(
    D = dist.mat,
    D1 = dist.mat.1,
    D2 = dist.mat.2,
    x1 = data.1$logPGA,
    x2 = data.2$logPGA,
    damage_ind_1 = damage_ind_1,
    damage_lower_1 = damage_lower_1,
    damage_state_1 = damage_state_1,
    damage_ind_2 = damage_ind_2,
    damage_lower_2 = damage_lower_2,
    damage_state_2 = damage_state_2,
    no_states = no_states,
    log_kappa = 0, # Fix smoothness parameter to INLA default to separate noise and smoothness.
    log_kappa1 = 0,
    log_kappa2 = 0
  )

  obj <- MakeADFun(
    data = input_data,
    parameters = parameters,
    random = c('field', 'field1', 'field2'),
    DLL = "spatialordinal_newmarginals_shared")

  its <- 150
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, trace = 0), upper = upper.lim, lower = lower.lim) #, abs.tol = 1e-3, rel.tol = 1e-2, x.tol = 1.5e-2, xf.tol = 2.2e-2)

  report <- obj$report()

  sd_out <- sdreport(obj, getJointPrecision = TRUE)

  sd.est <- diag(sd_out$cov.fixed)

  results <- opt
  results$field <- report$field
  results$field1 <- report$field1
  results$field2  <- report$field2
  results$par.sd <- sd.est

  return(results)

}
