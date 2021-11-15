#' Convert cut-off factors used in spatial ordinal model fit to values on the latent variable scale.
#'
#' @param par.est A vector with named components corresponding to the parameter estimates obtained from the \code{spatial_ordinal} function.
#' @return A vector with cut-off factors converted into cut-off values on the latent variable scale.
#' @export
#'
#' @examples
#' par_est <- c(1.8039257, -1.4570681, -16.6211310, -1.6039472,
#'              -0.9755298, -2.6170334, -0.8107559, -1.1138335,
#'              -1.4969662, -0.5405564, -0.8678417, 0.6900355,
#'              0.4677500, 0.3603201, 0.3383508, 0.3626860,
#'              0.4735808, 0.2284372, 0.4954409, 0.3912975,
#'              0.3595872, 0.3817603, 0.4913464)
#' names(par_est) <- c("log_phi", "log_sigma_2", "log_tau_2", "log_phi1",
#'                     "log_tau1_2", "log_sigma1_2", "log_phi2", "log_tau2_2",
#'                     "log_sigma2_2", "log_slope1", "log_slope2", "c_factor1",
#'                     "c_factor1", "c_factor1", "c_factor1", "c_factor1",
#'                     "c_factor1", "c_factor2", "c_factor2", "c_factor2",
#'                     "c_factor2", "c_factor2", "c_factor2")
#' convert_cutoffs(par_est)

convert_cutoffs <- function(par.est){

  cutoffs_factor1 <- par.est[names(par.est) == "c_factor1"]
  cutoffs1 <- cutoffs_factor1
  names(cutoffs1) <- rep("cutoffs1", length(cutoffs1))
  cutoffs_factor2 <- par.est[names(par.est) == "c_factor2"]
  cutoffs2 <- cutoffs_factor2
  names(cutoffs2) <- rep("cutoffs2", length(cutoffs2))

  for (i in 2:length(cutoffs1)){
    cutoffs1[i] <-  cutoffs1[i] + sum(cutoffs_factor1[1:(i-1)])
    cutoffs2[i] <-  cutoffs2[i] + sum(cutoffs_factor2[1:(i-1)])
  }

  same.par <- par.est[!(names(par.est) %in% c("c_factor1", "c_factor2"))]
  temp.names <- c(names(same.par), rep("cutoffs1", length(cutoffs1)), rep("cutoffs2", length(cutoffs2)))
  new.par <- c(same.par, cutoffs1, cutoffs2)
  names(new.par) <- temp.names
  return(new.par)
}
