#' Toy Multiplexed Imaging Data
#'
#' A small synthetic dataset mimicking multiplexed tissue imaging data, used in
#' package examples and tests.  The list contains one sample named
#' \code{"ToySample"}.
#'
#' @format A named list with one element:
#' \describe{
#'   \item{\code{ToySample}}{A \code{data.frame} with the following columns:
#'     \describe{
#'       \item{\code{x}}{Numeric. X coordinate in microns.}
#'       \item{\code{y}}{Numeric. Y coordinate in microns.}
#'       \item{\code{phenotype}}{Character. Cell phenotype label.  Values are
#'         \code{"B cell"}, \code{"T cell"}, and \code{"Other"}.}
#'     }
#'   }
#' }
#'
#' @source Synthetically generated for package examples.
#' @references Amiryousefi et al. (2025) \doi{10.1101/2025.09.21.677465}
#' @examples
#' data(toy_ldata)
#' str(toy_ldata)
#' table(toy_ldata[["ToySample"]]$phenotype)
"toy_ldata"
