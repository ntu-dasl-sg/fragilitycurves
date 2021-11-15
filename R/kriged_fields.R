#' Krige field estimates obtained via a spatial ordinal model fit.
#'
#' @param model.fit The output of a \code{spatial_ordinal} model fit.
#' @param data.1 The damage dataset of Building Category 1 used in the fit.
#' @param data.2 The damage dataset of Building Category 2 used in the fit.
#' @param zoom.shp A shapefile of the area to zoom into and bound the kriging by (optional).
#' @param PGA.raster A raster of the peak ground acceleration (PGA) which covers \code{zoom.shp} if applicable and whose resolution determines that of the outputs'.
#' @return A list of six rasters corresponding to the kriged shared field and building category specific fields together with their kriging variances.
#' @importFrom geoR krige.conv
#' @importFrom geoR krige.control
#' @import raster
#' @export
#'
#' @examples
#' library(raster)
#' library(rgdal)
#' library(sp)
#'
#' data(damage_simulation)
#' data(mean_PGA)
#' data(demo_spatial_fit)
#' data(haiti_admin2)
#'
#' data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
#' data.subset.2 <- damage_simulation[damage_simulation$building_cat == 2, ]
#'
#' # Obtain the shapefile for Port-au-Prince:
#'
#' pap_shp <- haiti_admin2[haiti_admin2$ADM2_EN == "Port-au-Prince", ]
#' pap_shp <- sp::spTransform(pap_shp, CRS("+proj=utm +zone=18 ellps=WGS84"))
#'
#' kriged_rasters <- kriged_fields(demo_spatial_fit, data.subset.1,
#'                                 data.subset.2, pap_shp, mean_PGA)


kriged_fields <- function(model.fit, data.1, data.2, zoom.shp = NULL, PGA.raster){


  data <- rbind(data.1, data.2)

  data$field <- model.fit$field; data.1$field1 <- model.fit$field1; data.2$field2 <- model.fit$field2

  data.1$field <- model.fit$field[1:nrow(data.1)]
  data.2$field <- model.fit$field[nrow(data.1)+1:nrow(data.2)]

  kappa.val <- kappa1.val <- kappa2.val <- 1

  if(!is.null(zoom.shp)){

    zoom.extent <- extent(zoom.shp)
    zoom.extent@xmin <- zoom.extent@xmin - 10000
    zoom.extent@ymin <- zoom.extent@ymin - 10000
    zoom.extent@xmax <- zoom.extent@xmax + 10000
    zoom.extent@ymax <- zoom.extent@ymax + 10000

    PGA.raster <- crop(PGA.raster, zoom.extent)

  }

  PGA.raster[is.na(PGA.raster)] <- -9999
  raster_pts <- rasterToPoints(PGA.raster, spatial = TRUE)
  coords_full <- raster_pts@coords/1000

  temp.time <- proc.time()[3]
  kriging_res <- geoR::krige.conv(coords = data[, c("Easting", "Northing")]/1000, data = data$field, locations = coords_full,
                                  krige = krige.control(type.krige = "sk", trend.d = "cte", trend.l = "cte",
                                                        obj.model = NULL, beta=0, cov.pars= as.numeric(c(exp(model.fit$par[names(model.fit$par) == "log_sigma_2"]), exp(model.fit$par[names(model.fit$par) == "log_phi"]))), kappa = kappa.val,
                                                        nugget = as.numeric(exp(model.fit$par[names(model.fit$par) == "log_tau_2"])), micro.scale = 0, dist.epsilon = 1e-20, lambda = 1))
  time.taken.3 <- proc.time()[3] - temp.time # 29.02s.

  res_df <- cbind(raster_pts@coords, kriging_res$predict)
  var_df <- cbind(raster_pts@coords, kriging_res$krige.var)
  field_ras <- rasterFromXYZ(res_df)
  var_ras <-  rasterFromXYZ(var_df)


  temp.time <- proc.time()[3]
  kriging_res_1 <- krige.conv(coords = data.1[, c("Easting", "Northing")]/1000, data = data.1$field1, locations = coords_full,
                              krige = krige.control(type.krige = "sk", trend.d = "cte", trend.l = "cte",
                                                    obj.model = NULL, beta=0, cov.pars= as.numeric(c(exp(model.fit$par[names(model.fit$par) == "log_sigma1_2"]), exp(model.fit$par[names(model.fit$par) == "log_phi1"]))), kappa = kappa1.val,
                                                    nugget = as.numeric(exp(model.fit$par[names(model.fit$par) == "log_tau_2"])), micro.scale = 0, dist.epsilon = 1e-10, lambda = 1))
  time.taken.4 <- proc.time()[3] - temp.time # 5.73s.

  temp.time <- proc.time()[3]
  kriging_res_2 <- krige.conv(coords = data.2[, c("Easting", "Northing")]/1000, data = data.2$field2, locations = coords_full,
                              krige = krige.control(type.krige = "sk", trend.d = "cte", trend.l = "cte",
                                                    obj.model = NULL, beta=0, cov.pars= as.numeric(c(exp(model.fit$par[names(model.fit$par) == "log_sigma2_2"]), exp(model.fit$par[names(model.fit$par) == "log_phi2"]))), kappa = kappa2.val,
                                                    nugget = as.numeric(exp(model.fit$par[names(model.fit$par) == "log_tau_2"])), micro.scale = 0, dist.epsilon = 1e-10, lambda = 1))
  time.taken.5 <- proc.time()[3] - temp.time # 13.22s.

  res1_df <- cbind(raster_pts@coords, kriging_res_1$predict)
  var1_df <- cbind(raster_pts@coords, kriging_res_1$krige.var)
  field1_ras <- rasterFromXYZ(res1_df)
  var1_ras <-  rasterFromXYZ(var1_df)

  res2_df <- cbind(raster_pts@coords, kriging_res_2$predict)
  var2_df <- cbind(raster_pts@coords, kriging_res_2$krige.var)
  field2_ras <- rasterFromXYZ(res2_df)
  var2_ras <-  rasterFromXYZ(var2_df)

  return(list("field_ras" = field_ras, "field1_ras" = field1_ras, "field2_ras" = field2_ras, "var_ras" = var_ras, "var1_ras" = var1_ras, "var2_ras" = var2_ras))

}
