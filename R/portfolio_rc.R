#' Compute mean replacement cost for building portfolio per event based on given cut-off values and the standard deviation of an ordinal model.
#'
#' @param cutoffs A vector of cutoff values (sorted in increasing order) for the ordinal model.
#' @param damage.states A vector of the damage state names e.g. their central damage factor values (character or factor).
#' @param lv A matrix containing latent variable means per event (by column) and per building in portfolio (by row).
#' @param ordinal.sd Standard deviation of the ordinal model. This is the denominator in the standard normal cumulative distribution function used to calculate exceedance probabilities.
#' @param replacement.cost The replacement cost of a building in the portfolio.
#' @param no.building A vector of the number of buildings corresponding to each row of \code{lv} if its rows correspond to grid locations instead of individual buildings.
#' @return A vector containing mean replacement costs per event for the building portfolio.
#' @importFrom stats pnorm
#' @export
#'
#' @examples
#' data(demo_spatial_fit)
#'
#' # Convert the cut-off factors to cut-off values on the latent variable scale:
#' new_par <- convert_cutoffs(demo_spatial_fit$par)
#'
#' # Spatial model parameters:
#' field.sigma2 <- exp(new_par["log_sigma_2"]);
#' field.tau2 <- exp(new_par["log_tau_2"]);
#' field1.sigma2 <- exp(new_par["log_sigma1_2"]);
#' field1.tau2 <- exp(new_par["log_tau1_2"]);
#' beta1 <- exp(new_par["log_slope1"]);
#' beta2 <- exp(new_par["log_slope2"]);
#'
#' data(oq_sim)
#'
#' # The following takes 2 minutes on a PC with characteristics:
#' # Intel(R) Xeon (R) W-2112 CPU Processor @ 3.60GHz; 32GB of RAM; Windows 10 64-bit:
#' nonspat_lv <- lv_sim(model = "Non-spatial", data = oq_sim,
#'                      fieldsim = NULL, field1sim = NULL, field2sim = NULL,
#'                      slope1 = beta1, slope2 = beta2,
#'                      shared_sill = field.tau2 + field.sigma2)
#'
#' cutoffs1 <- new_par[names(new_par) == "cutoffs1"];
#' CDF_breaks <- sort(unique(data.subset.1$CDF), decreasing = FALSE)
#'
#' # The following takes 4 minutes on a PC with characteristics:
#' # Intel(R) Xeon (R) W-2112 CPU Processor @ 3.60GHz; 32GB of RAM;
#' # nonspat.rc1 <- portfolio_rc(cutoffs1, CDF_breaks, nonspat_lv$lv1,
#' #                             sqrt(field1.tau2 + field1.sigma2),
#' #                             replacement.cost = 1)

portfolio_rc <- function(cutoffs, damage.states, lv, ordinal.sd, replacement.cost = 1, no.building = NULL){
  if(is.null(no.building)){no.building <- rep(1, nrow(lv))}
  m.rc <- matrix(NA, nrow = nrow(lv), ncol = ncol(lv))
  no_states <- length(cutoffs) + 1
  for (i in 1:ncol(lv)){ # For each event.
    temp <- matrix(0, nrow = nrow(lv), ncol = no_states) # Compute the damage probabilities per state.
    lv_i <- lv[, i]
    for (j in 1:(no_states-2)){ # Compute the probability of each damage state for all locations.
      temp[, j+1] <-   pnorm((cutoffs[j+1] - lv_i)/ordinal.sd) - pnorm((cutoffs[j] - lv_i)/ordinal.sd)
    }
    temp[, no_states] <- 1 - pnorm((cutoffs[no_states-1] - lv_i)/ordinal.sd)
    temp[, 1] <- 1 - rowSums(temp)
    # Compute the mean replacement cost as a weighted average:
    m.rc[, i] <- apply(temp, MARGIN = 1, FUN = function(x){sum(x*as.numeric(as.character(damage.states)))/100})
  }
  rc <- colSums(matrix(rep(no.building), ncol(lv), nrow = nrow(lv), ncol = ncol(lv))*m.rc)*replacement.cost
  return(rc)
}
