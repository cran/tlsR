#' tlsR: Detection and Spatial Analysis of Tertiary Lymphoid Structures
#'
#' Fast, reproducible detection and quantitative analysis of tertiary lymphoid
#' structures (TLS) in multiplexed tissue imaging data.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Load or prepare a named list of data frames (\code{ldata}), one per
#'         tissue sample. Each data frame must contain columns \code{x}, \code{y}
#'         (spatial coordinates in microns), and \code{phenotype} (character:
#'         \code{"B cell"} / \code{"T cell"} / other).
#'   \item Run \code{\link{detect_TLS}} to label B+T co-localised regions.
#'   \item (Optional) Run \code{\link{scan_clustering}} to identify windows of
#'         significant immune clustering via local Ripley's L.
#'   \item Run \code{\link{calc_icat}} to score the internal linearity/organisation
#'         of each detected TLS.
#'   \item Run \code{\link{detect_tic}} to identify T-cell clusters outside TLS.
#'   \item Use \code{\link{summarize_TLS}} to obtain a tidy summary table.
#'   \item Use \code{\link{plot_TLS}} to produce publication-ready spatial plots.
#' }
#'
#' @references
#' Amiryousefi et al. (2025) \doi{10.1101/2025.09.21.677465}
#'
#' @docType package
#' @name tlsR-package
#' @aliases tlsR
#' @importFrom methods is
"_PACKAGE"
