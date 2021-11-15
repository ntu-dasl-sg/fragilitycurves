#' Compute and plot the exceedance probability rasters for the different damage states.
#'
#' @param category Set this to 1 (Building category 1) or 2 (Building category 2).
#' @param new.par The vector with parameter estimates corresponding to a \code{spatial_ordinal} model fit or its equivalent with cut-off factors converted to the latent variable scale via \code{convert_cutoffs}.
#' @param damage.states A vector of the damage state names e.g. their central damage factor values (character or factor).
#' @param latent.raster A raster of the latent variable mean surface. This could be the output of the \code{latent_var} function.
#' @param study.shp A shapefile of the study area which covers the area spanned by \code{PGA.raster} (optional). If not provided, the plots of the exceedance probabilities will not be computed.
#' @return A dataframe with columns corresponding to "Easting", "Northing" and the exceedance probabilities of the damage states.
#' @importFrom gridExtra grid.arrange
#' @importFrom stats pnorm
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
#' study_shp <- sp::spTransform(haiti_admin0, CRS("+proj=utm +zone=18 ellps=WGS84"))
#'
#' # Obtain the shapefile for Port-au-Prince:
#' pap_shp <- haiti_admin2[haiti_admin2 == "Port-au-Prince", ]
#' pap_shp <- sp::spTransform(pap_shp, CRS("+proj=utm +zone=18 ellps=WGS84"))
#'
#' data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
#'
#' kriged_rasters <- kriged_fields(demo_spatial_fit, data.subset.1,
#'                                 data.subset.2, pap_shp, mean_PGA)
#' latent_var_1 <- latent_var(category = 1, new_par, kriged_rasters,
#'                            mean_PGA, study_shp)
#'
#' CDF_breaks <- sort(unique(data.subset.1$CDF), decreasing = FALSE)
#'
#' exceed_prob_1 <- prob_exceed(1, new_par, CDF_breaks, latent_var_1, study_shp)

prob_exceed <- function(category, new.par, damage.states, latent.raster, study.shp){

  if(category == 1){
    cutoffs <- new.par[names(new.par) == "cutoffs1"]
    marginal.sd <- sqrt(exp(new.par[names(new.par) == "log_tau1_2"]))
  }else{
    cutoffs <- new.par[names(new.par) == "cutoffs2"]
    marginal.sd <- sqrt(exp(new.par[names(new.par) == "log_tau2_2"]))
  }

  raster.extent <- latent.raster@extent

  if (!is.null(study.shp)){
    p.list <- list()
    study.shp_crop <- crop(study.shp, raster.extent)
    study.reg.map <- geom_polygon(data = study.shp_crop, aes(x = long, y = lat, group = group), colour = "black", fill = NA)
  }


  for (i in 1:length(cutoffs)){

    ex_prob <- 1 - pnorm((cutoffs[i] - raster::as.matrix(latent.raster))/marginal.sd)
    ex_prob_raster <- raster::raster(ex_prob, xmn = raster.extent[1], xmx = raster.extent[2], ymn = raster.extent[3], ymx = raster.extent[4],crs = crs(latent.raster))

    ex.raster.df <- as.data.frame(raster::rasterToPoints(ex_prob_raster))

    colnames(ex.raster.df) <- c("x", "y", "ex.prob")

    if (!is.null(study.shp)){
      p.list[[i]] <- ggplot() + geom_raster(data = ex.raster.df, aes(x=x, y=y, fill = ex.prob)) + coord_equal() + labs(x = "", y = "", fill = "P(Exceedance)") + ggtitle(paste("State ", i, sep = "")) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + study.reg.map + scale_fill_viridis_c(limits = c(0, 1))
    }

    if(i == 1){ex.prob.df <- ex.raster.df}else{ex.prob.df[, as.character(damage.states[i])] <- ex.raster.df$ex.prob}

  }

  if (!is.null(study.shp)){
    do.call(grid.arrange, c(p.list, nrow = 2))
  }
  names(ex.prob.df) <- c("Easting", "Northing", as.character(damage.states[-length(damage.states)]))

  return(ex.prob.df)

}
