#' Excerpt from the Open Exoplanet Catalogue Data Set
#'
#' Data fields include planet and host star attributes.
#'
#' The variables are as follows:
#'
#' \itemize{
#' \item semimajoraxis: Semi-major axis of the planetary object's orbit in astronomical units
#' \item mass: Mass of the planetary object in Jupiter masses
#' \item radius: Radius of the planetary object in Jupiter radii
#' \item period: Orbital period of the planetary object in days
#' \item eccentricity: Eccentricity of the planetary object's orbit
#' \item hoststar_mass: Mass of the host star in solar masses
#' \item hoststar_radius: Radius of the host star in solar radii
#' \item hoststar_metallicity: Metallicity of the host star
#' \item hoststar_temperature: Effective temperature of the host star in Kelvin
#' \item binaryflag: Flag indicating the type of planetary system
#' }
#'
#' @docType data
#' @keywords datasets
#' @name exoplanet
#' @usage data(exoplanet)
#' @format A data frame with 223 rows and 11 variables
#' @source Dataset downloaded from the Open Exoplanet Catalogue Repository.
#' \url{https://github.com/OpenExoplanetCatalogue/oec_tables/}
#'
#' Creators:
#'
#' 1. Prof. Hanno Rein, Department for Physical and Environmental Sciences.
#' University of Toronto at Scarborough
#' Toronto, Ontario M1C 1A4
#' hanno.rein 'at' utoronto.ca
#'
"exoplanet"
