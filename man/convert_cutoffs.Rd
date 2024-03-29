% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_cutoffs.R
\name{convert_cutoffs}
\alias{convert_cutoffs}
\title{Convert cut-off factors used in spatial ordinal model fit to values on the latent variable scale.}
\usage{
convert_cutoffs(par.est)
}
\arguments{
\item{par.est}{A vector with named components corresponding to the parameter estimates obtained from the \code{spatial_ordinal} function.}
}
\value{
A vector with cut-off factors converted into cut-off values on the latent variable scale.
}
\description{
Convert cut-off factors used in spatial ordinal model fit to values on the latent variable scale.
}
\examples{
par_est <- c(1.8039257, -1.4570681, -16.6211310, -1.6039472,
             -0.9755298, -2.6170334, -0.8107559, -1.1138335,
             -1.4969662, -0.5405564, -0.8678417, 0.6900355,
             0.4677500, 0.3603201, 0.3383508, 0.3626860,
             0.4735808, 0.2284372, 0.4954409, 0.3912975,
             0.3595872, 0.3817603, 0.4913464)
names(par_est) <- c("log_phi", "log_sigma_2", "log_tau_2", "log_phi1",
                    "log_tau1_2", "log_sigma1_2", "log_phi2", "log_tau2_2",
                    "log_sigma2_2", "log_slope1", "log_slope2", "c_factor1",
                    "c_factor1", "c_factor1", "c_factor1", "c_factor1",
                    "c_factor1", "c_factor2", "c_factor2", "c_factor2",
                    "c_factor2", "c_factor2", "c_factor2")
convert_cutoffs(par_est)
}
