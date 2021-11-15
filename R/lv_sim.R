#' Compute latent variable means for the two building categories for the spatial ordinal model or its submodels.
#'
#' @param model This takes the values "Non-spatial" for the non-spatial submodel, "IM-spatial" for the submodel with the shared spatial field but no building category specific field, and "Damage-spatial" for the full spatial ordinal model.
#' @param data A dataframe with columns "lon", "lat" and "building cat" denoting the longitude, latitude and building category ("1"/"2") of the buildings. Appended to the right of this are columns with event IDs as names, containing the simulated log(PGA) values. 'NA' values are treated as -Inf.
#' @param fieldsim A dataframe with simulated shared field values for the buildings and events in \code{data}. This is required unless \code{model = "Non-spatial"}.
#' @param field1sim A dataframe with simulated building category 1 field values for the category 1 buildings and events in \code{data}. This is only required for \code{model = "Damage-spatial"}.
#' @param field2sim A dataframe with simulated building category 2 field values for the category 2 buildings and events in \code{data}. This is only required for \code{model = "Damage-spatial"}.
#' @param slope1 Slope coefficient for the log(PGA) component of the latent variable mean formula for buildling category 1.
#' @param slope2  Slope coefficient for the log(PGA) component of the latent variable mean formula for buildling category 2.
#' @param shared_sill Sill of the shared spatial field. This is required unless \code{model = "Damage-spatial"}.
#' @return A list of two matrices containing latent variable means per event (by column) and per building (by row): \code{lv1} for Building category 1 and \code{lv2} for Building category 2.
#' @importFrom stats rnorm
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
#' beta1 <- exp(new_par["log_slope1"]);
#' beta2 <- exp(new_par["log_slope2"]);
#'
#' data(oq_sim)
#'
#' # The following takes 2 minutes on a PC with characteristics:
#' # Intel(R) Xeon (R) W-2112 CPU Processor @ 3.60GHz; 32GB of RAM; Windows 10 64-bit:
#' #nonspat_lv <- lv_sim(model = "Non-spatial", data = oq_sim,
#' #                    fieldsim = NULL, field1sim = NULL, field2sim = NULL,
#' #                    slope1 = beta1, slope2 = beta2,
#' #                    shared_sill = field.tau2 + field.sigma2)

lv_sim <- function(model = "Non-spatial", data = subset_data, fieldsim = NULL, field1sim = NULL, field2sim = NULL, slope1 = beta1, slope2 = beta2, shared_sill = field.tau2 + field.sigma2){

  logPGAsim <- data[, !(colnames(data) %in% c("lon", "lat", "building_cat"))]
  logPGAsim[is.na(logPGAsim)] <- -Inf

  nonspat.lv1 <- as.matrix(slope1*logPGAsim)[data$building_cat == 1, ]
  nonspat.lv2 <- as.matrix(slope2*logPGAsim)[data$building_cat == 2, ]

  if(model == "IM-spatial" | model == "Damage-spatial"){

    sharedspat.lv1 <- nonspat.lv1 + as.matrix(slope1*fieldsim@data)[data$building_cat == 1, ]
    sharedspat.lv2 <- nonspat.lv2 + as.matrix(slope2*fieldsim@data)[data$building_cat == 2, ]

    if(model == "Damage-spatial"){
      damagespat.lv1 <- sharedspat.lv1 + as.matrix(field1sim@data)
      damagespat.lv2 <- sharedspat.lv2 + as.matrix(field2sim@data)
    }
  }

  if(model == "Non-spatial"){

    # Generate random error terms to be added to the non-spatial latent variable mean only.
    set.seed(5)
    iid.error3 <- matrix(rnorm(n = nrow(data)*ncol(logPGAsim), mean = 0, sd = sqrt(shared_sill)), nrow = nrow(data), ncol = ncol(logPGAsim))

    lv1 <- nonspat.lv1 + slope1*iid.error3[data$building_cat == 1, ]
    lv2 <- nonspat.lv2 + slope2*iid.error3[data$building_cat == 2, ]
  }
  if(model == "IM-spatial"){
    lv1 <- sharedspat.lv1
    lv2 <- sharedspat.lv2
  }
  if(model == "Damage-spatial"){
    lv1 <- damagespat.lv1
    lv2 <- damagespat.lv2
  }

  return(list("lv1" = lv1, "lv2" = lv2))
}
