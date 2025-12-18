#' Toy multiplexed imaging data for examples
#'
#' A small synthetic dataset mimicking multiplexed tissue imaging data.
#' It contains one sample named "ToySample" with columns required by tlsR functions.
#'
#' @format A named list with one element:
#' \describe{
#'   \item{ToySample}{A data frame with columns:
#'     \itemize{
#'       \item x: x-coordinate in microns
#'       \item y: y-coordinate in microns
#'       \item coarse_phen_vec: Cell phenotype ("B cells", "T cells", or "Other")
#'       \item row_index: Integer row index (1 to nrow)
#'       \item cflag: Integer flag column (0 for all cells)
#'     }
#'   }
#' }
#' @examples
#' data(toy_ldata)
#' str(toy_ldata[["ToySample"]])
#' plot(toy_ldata[["ToySample"]]$x, toy_ldata[["ToySample"]]$y,
#'      col = as.factor(toy_ldata[["ToySample"]]$coarse_phen_vec),
#'      pch = 19, cex = 0.5, main = "Toy sample cells")
"toy_ldata"
