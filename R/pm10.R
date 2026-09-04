#' Realistic data which illustrate the usage of the package Stem
#'
#' @name pm10
#' @aliases pm10
#' @docType data
#'
#' @description
#' This simple dataset is a list of three objects and refers to 22 spatial locations and 366 time points.
#'
#' @usage data(pm10)
#'
#' @format A list with three objects containing the following components:
#' \describe{
#'   \item{coords}{A matrix with the coordinates of the 22 spatial locations.}
#'   \item{covariates}{An 8052 by 3 matrix referring to the following covariates: \emph{intercept, emissions (g/s), and altitude (km)}. The first 366 rows refer to the first spatial
#'   location, the rows from 367 to 732 refer to the second spatial location, and so on.}
#'   \item{z}{A 366 by 22 observation matrix referring to \emph{PM10 concentration measurements} (log scale).}
#' }
#'
#' @references
#' Fasso, A., Cameletti, M. (2007) \emph{A general spatio-temporal model for environmental data}. Tech.rep. n.27 \emph{Graspa} - The Italian Group of Environmental Statistics.
#'
#' @author Michela Cameletti \email{michela.cameletti@unibg.it}
#'
#' @examples
#' data(pm10)
#' names(pm10)
#'
#' # Plot the coordinates
#' dim(pm10$coords)
#' plot(pm10$coords[,1], pm10$coords[,2], xlab=colnames(pm10$coords)[1],
#'      ylab=colnames(pm10$coords)[2])
#'
#' # Plot the data
#' dim(pm10$z)
#'
#' # Summary by station
#' apply(pm10$z, 2, summary)
#'
#' # Plot the time series for station n.22
#' plot(pm10$z[,22], type="l", xlab="Days", ylab="PM10 concentrations (log)")
#'
#' # Plot the station altitude
#' plot(pm10$covariates[,3], ylab=colnames(pm10$covariates)[3], xaxt="n", xlab="")
#' positions <- seq(1, 8052, 366) + 366/2
#' axis(1, at=positions, labels=rownames(pm10$coords), las=2)
#'
#' @keywords datasets
"pm10"
