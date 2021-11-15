#' Compute the latent variable mean raster for a building category and plot its contributing terms.
#'
#' @param category Set this to 1 (Building category 1) or 2 (Building category 2).
#' @param new.par The vector with parameter estimates corresponding to a \code{spatial_ordinal} model fit or its equivalent with cut-off factors converted to the latent variable scale via \code{convert_cutoffs}.
#' @param field.rasters A list of rasters where the first three correspond to the kriged shared field and two building category specific fields respectively. This could be the output of \code{kriged_fields}.
#' @param PGA.raster A raster of the peak ground acceleration (PGA) which covers the extent of the rasters in \code{field.rasters}.
#' @param study.shp A shapefile of the study area which covers the area spanned by \code{PGA.raster} (optional). If not provided, the plot of the contributing terms will not be computed.
#' @return A raster of the latent variable mean surface computed according to the fitted spatial ordinal model for the building category.
#' @importFrom gridExtra grid.arrange
#' @import raster
#' @import ggplot2
#' @export
#'
#' @examples
#' library(raster)
#' library(rgdal)
#' library(sp)
#' library(rgeos)
#'
#' data(damage_simulation)
#' data(mean_PGA)
#' data(demo_spatial_fit)
#' data(haiti_admin2)
#' data(haiti_admin0)
#'
#' # Convert the cut-off factors to cut-off values on the latent variable scale:
#' new_par <- convert_cutoffs(demo_spatial_fit$par)
#'
#' study_shp <- sp::spTransform(haiti_admin0, sp::CRS("+proj=utm +zone=18 ellps=WGS84"))
#'
#' # Obtain the shapefile for Port-au-Prince:
#' pap_shp <- haiti_admin2[haiti_admin2 == "Port-au-Prince", ]
#' pap_shp <- sp::spTransform(pap_shp, sp::CRS("+proj=utm +zone=18 ellps=WGS84"))
#'
#' data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
#'
#' kriged_rasters <- kriged_fields(demo_spatial_fit, data.subset.1,
#'                                 data.subset.2, pap_shp, mean_PGA)
#' latent_var_1 <- latent_var(category = 1, new_par, kriged_rasters,
#'                            mean_PGA, study_shp)

latent_var <- function(category, new.par, field.rasters, PGA.raster, study.shp = NULL){


  field_ras <- field.rasters[[1]];
  if(category == 1){
    cat_field_ras <- field.rasters[[2]]
    slope <- exp(new.par[names(new.par) == "log_slope1"])
  }else{
    cat_field_ras <- field.rasters[[3]]
    slope <- exp(new.par[names(new.par) == "log_slope2"])
  }

  raster.extent <- field_ras@extent

  PGA.raster <- raster::crop(PGA.raster, raster.extent)

  term1 <- slope * log(PGA.raster)
  term2 <- slope *  field_ras

  latent.var <- term1 + term2 + cat_field_ras

  if(!is.null(study.shp)){

    # Use ggplot2 for the images:

    raster.3.df <- as.data.frame(raster::rasterToPoints(cat_field_ras))
    colnames(raster.3.df) <- c("x", "y", "cat.field")

    study.shp_crop <- raster::crop(study.shp, raster.extent)
    study.reg.map <- geom_polygon(data = study.shp_crop, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

    raster_1 <- term1
    raster.1.df <- as.data.frame(raster::rasterToPoints(raster_1))
    colnames(raster.1.df) <- c("x", "y", "beta.lPGA")

    lPGA.map <- ggplot(data = raster.1.df, aes(x=x, y=y)) + geom_raster(aes(fill = beta.lPGA)) + scale_fill_viridis_c()  + coord_equal() + labs(x = "Easting", y = "Northing", fill = "")
    map.1 <- lPGA.map + study.reg.map + ggtitle(expression(paste("(a) ", beta, "log(PGA)"))) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

    raster_2 <- term2
    raster.2.df <- as.data.frame(rasterToPoints(raster_2))
    colnames(raster.2.df) <- c("x", "y", "beta.field")

    field.map <- ggplot(data = raster.2.df, aes(x=x, y=y)) + geom_raster(aes(fill = beta.field)) + scale_fill_viridis_c() + coord_equal() + labs(x = "Easting", y = "", fill = "")
    map.2 <- field.map + study.reg.map + ggtitle(expression(paste("(b) ", beta*u^{S}))) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

    cat.field.map <- ggplot(data = raster.3.df, aes(x=x, y=y)) + geom_raster(aes(fill = cat.field)) + scale_fill_viridis_c() + coord_equal() + labs(x = "Easting", y = "", fill = "")
    map.3 <- cat.field.map + study.reg.map + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

    if(category == 1){map.3 <- map.3 + ggtitle(expression(paste("(c) ", u^{(1)})))}else{map.3 <- map.3 + ggtitle(expression(paste("(c) ", u^{(2)})))}

    gridExtra::grid.arrange(map.1, map.2, map.3, ncol = 3)

  }

  return(latent.var)

}
