#' Plot the fitted variograms for the shared and building category specific spatial fields.
#'
#' @param par.est A vector with named components corresponding to the parameter estimates obtained from the \code{spatial_ordinal} function.
#' @param shared.range A vector containing the values at which the variogram of the shared spatial field should be evaluated.
#' @param cat.range A vector containing the values at which the variograms of the building category specific spatial fields should be evaluated.
#' @return A plot of three variograms corresponding to the shared and the two building category specific spatial fields.
#' @importFrom geoR matern
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
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
#' shared_range <- seq(0, 50, by = 0.2); cat_range <- seq(0, 5, by = 0.02)
#' vgm_plot(par_est, shared_range, cat_range)

vgm_plot <- function(par.est, shared.range, cat.range){

  shared.vgm <- as.numeric(exp(par.est["log_sigma_2"]))*(1-geoR::matern(shared.range, phi= as.numeric(exp(par.est["log_phi"])), kappa = 1)) + as.numeric(exp(par.est["log_tau_2"]))
  cat1.vgm <- as.numeric(exp(par.est["log_sigma1_2"]))*(1-geoR::matern(cat.range, phi= as.numeric(exp(par.est["log_phi1"])) , kappa = 1)) + as.numeric(exp(par.est["log_tau1_2"])) + as.numeric(exp(par.est["log_tau_2"]))

  cat2.vgm <- as.numeric(exp(par.est["log_sigma2_2"]))*(1-geoR::matern(cat.range, phi= as.numeric(exp(par.est["log_phi2"])) , kappa = 1)) + as.numeric(exp(par.est["log_tau2_2"])) + as.numeric(exp(par.est["log_tau_2"]))

  shared.df <- data.frame("x" = shared.range, "y" = shared.vgm)
  cat1.df <- data.frame("x" = cat.range, "y" = cat1.vgm)
  cat2.df <- data.frame("x" = cat.range, "y" = cat2.vgm)

  shared.plot <- ggplot(data = shared.df, aes(x = x, y = y)) + geom_line() + labs(x = "km", y = expression(gamma(h))) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(a) Shared spatial field")

  cat1.plot <- ggplot(data = cat1.df, aes(x = x, y = y)) + geom_line() + labs(x = "km", y = expression(gamma(h))) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(b) Building category 1 specific")

  cat2.plot <- ggplot(data = cat2.df, aes(x = x, y = y)) + geom_line() + labs(x = "km", y = expression(gamma(h))) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(c) Building category 2 specific")

  gridExtra::grid.arrange(shared.plot, cat1.plot, cat2.plot, ncol = 1)

}
