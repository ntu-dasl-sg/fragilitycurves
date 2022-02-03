#' Simulated damage for an earthquake event in Haiti.
#'
#' A dataset containing a simulation of damage information for 350 buildings in Haiti (175 unreinforced block walls and 175 stone masonry buildings in Haiti with 25 buildings in each damage grade in the damage scale used in ATC-13).
#'
#' @format A dataframe with 350 rows and 6 variables:
#' \describe{
#' \item{building_cat}{Building category, 1 for unreinforced block walls and 2 for stone masonry}
#' \item{CDF}{Central damage factor corresponding to the damage grade: 0 for Grade 1, 0.5 for Grade 2, 5 for Grade 3, 20 for Grade 4, 45 for Grade 5, 80 for Grade 6, 100 for Grade 7}
#' \item{logPGA}{Log-transformed peak ground acceleration (PGA) where PGA is given in %g}
#' \item{PGA}{PGA in %g}
#' \item{Easting}{Easting coordinate of the building}
#' \item{Northing}{Northing coordinate of the building}
#' }
#' @source Simulated from a spatial ordinal model using logPGA simulations from OpenQuake.
"damage_simulation"

#' Spatial ordinal model fit to the simulated damage for an earthquake event in Haiti.
#'
#' A list containing the output of \code{spatial_ordinal} when the function is applied to the damage_simulation data.
#'
#' @format A list with 10 variables:
#' \describe{
#' \item{par}{A vector containing the estimates for for the model parameters: \code{c_factor1} (first cut-off value and cut-off ratios for Building Category 1), \code{c_factor2} (first cut-off value and cut-off ratios for Building Category 2), \code{log_phi} (log-transformed range parameter for the shared spatial field), \code{log_sigma_2} (log-transformed partial sill parameter for the shared spatial field), \code{log_phi1} (log-transformed range parameter for the Building Category 1 spatial field), \code{log_tau1_2} (log-transformed nugget parameter for the Building Category 1 spatial field), \code{log_sigma1_2} (log-transformed partial sill parameter for the Building Category 1 spatial field), \code{log_phi2} (log-transformed range parameter for the Building Category 2 spatial field), \code{log_tau2_2} (log-transformed nugget parameter for the Building Category 2 spatial field), \code{log_sigma2_2} (log-transformed partial sill parameter for the Building Category 2 spatial field), \code{log_slope1} (log-transformed slope parameter for Building Category 1) and \code{log_slope2} (log-transformed slope parameter for Building Category 2).}
#' \item{objective}{Final value of the objective function e.g. negative log-likelihood.}
#' \item{convergence}{Convergence code of the nlminb algorithm. 0 represents convergence.}
#' \item{iterations}{Number of iterations taken in optimisation procedure before stopping. Note that a maximum of 150 iterations have been set in the \code{spatial_ordinal} function.}
#' \item{evaluations}{Number of function and gradient evaluations taken in optimisation procedure before stopping.}
#' \item{message}{Message obtained from the nlminb algorithm.}
#' \item{field}{Estimated values of the field shared between different building categories with long-range spatial correlation.}
#' \item{field1}{Estimated values of the spatial field associated with building category 1 with shorter-range spatial correlation.}
#' \item{field2}{Estimated values of the spatial field associated with building category 2 with shorter-range spatial correlation.}
#' \item{par.sd}{Estimated standard deviation obtained from the \code{TMB::sdreport} function for the parameters in \code{par}.}
#' }
#' @source Obtained from a spatial ordinal model fit to the simulated damage for an earthquake event in Haiti.
"demo_spatial_fit"

#' Simulated peak ground acceleration in %g for Haiti.
#'
#' A raster of mean peak ground acceleration in %g for Haiti.
#'
#' @format A 212 x 212 raster of 1 degree longitude/latitude resolution projected to the UTM coordinate system.
#' @source Simulated from a bespoke ground shaking model.
"mean_PGA"

#' Shapefiles of Haiti level 2 administrative units.
#'
#' A SpatialPolygonsDataFrame containing shapefiles of level 2 administrative units in Haiti.
#'
#' @format SpatialPolygonsDataFrame containing 140 shapefiles of level 2 administrative units in Haiti.
#' @source Haiti Subnational Administrative Boundaries contributed by the UN office for the Coordination of Humanitarian Affairs (OCHA) country office in Haiti, available on the Humanitarian Data Exchange (HDX).
"haiti_admin2"

#' Shapefile for the Haiti level 0, i.e. national administrative unit.
#'
#' A SpatialPolygonsDataFrame containing shapefiles of the level 0, i.e. national administrative unit in Haiti.
#'
#' @format SpatialPolygonsDataFrame containing 1 shapefile of the level 0, i.e. national administrative unit in Haiti.
#' @source Haiti Subnational Administrative Boundaries contributed by the UN office for the Coordination of Humanitarian Affairs (OCHA) country office in Haiti, available on the Humanitarian Data Exchange (HDX).
"haiti_admin0"

#' Simulated log-transformed peak ground acceleration (PGA) values for 300 buildings in Port-au-Prince, Haiti
#'
#' A dataframe of simulated log-transformed peak ground acceleration (PGA) values for 150 unreinforced block walls (Building Category 1) and 150 stone masonry (Building Category 2) buildings within a 2km x 2km grid in Port-au-Prince, Haiti, based on 1 million stochastic event sets of one year duration.
#'
#' @format A dataframe with 300 rows and 156494 variables:
#' \describe{
#' \item{lon}{Longitude coordinates (in degrees) of building.}
#' \item{lat}{Latitude coordinates (in degrees) of building.}
#' \item{building_cat}{Building category, 1 for unreinforced block walls and 2 for stone masonry.}
#' \item{`10`}{Simulated log(PGA) values for each building from simulated event 10. The remaining columns denote the simulated values from events where shaking was observed.}
#' }
#' @source Simulated using the OpenQuake Engine.
"oq_sim"

#' Matching list for simulated OpenQuake events and their stochastic event sets of one year duration.
#'
#' A dataframe of OpenQuake event IDs and their corresponding realisation ID denoting the year of occurrence or stochastic event set.
#'
#' @format A dataframe with 166882 rows and 2 variables:
#' \describe{
#' \item{event_id}{Event ID of a simulated earthquake event affecting the selected portfolio of buildings in \code{oq_sim} where the corresponding logPGA values are given in the column named after the event ID.}
#' \item{rlz_id}{Realisation ID indicating the year of occurrence of stochastic event set (SES) of one year duration. Note that we've considered a total of 1 million SES in our simulations.}
#' }
#' @source Simulated using the OpenQuake Engine.
"event_ses_list"
