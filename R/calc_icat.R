#' Calculate ICAT (Immune Cell Arrangement Trace) Index
#'
#' Quantifies the spatial spread and linear organisation of cells within a
#' detected TLS.  FastICA is applied to the (x, y) coordinates of TLS cells
#' to estimate independent components; the mixing matrix is used to
#' reconstruct the data, and the ICAT index is defined as the normalised
#' trace-standard-deviation of the reconstructed coordinates.
#'
#' @details
#' The ICAT index is computed as follows:
#' \enumerate{
#'   \item Centre the (x, y) coordinates of TLS cells.
#'   \item Run \code{fastICA} with 2 components.
#'   \item Reconstruct \eqn{\hat{X} = S A^T + \mu}.
#'   \item Let \eqn{v_1, v_2} be the marginal variances of \eqn{\hat{X}}.
#'   \item \deqn{\text{ICAT} = 100 \times
#'         \frac{\sqrt{v_1 + v_2 + 2\sqrt{v_1 v_2}}}{\text{nrow}(X)}}
#' }
#' This formulation is always non-negative: it measures the average
#' spatial spread (in microns) per cell rather than the signed trace of the
#' mixing matrix, which can be negative due to ICA sign ambiguity.
#' Higher ICAT values indicate a more spatially extended (structured) cluster.
#'
#' If the requested TLS contains fewer than 3 cells, or FastICA does not
#' converge, the function returns \code{NA_real_} with an informative message
#' rather than throwing an error.
#'
#' @param patientID Character. Sample name in \code{ldata}.
#' @param tlsID Numeric or integer. TLS identifier (value of
#'   \code{tls_id_knn} for the TLS of interest).
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#'
#' @return A single non-negative numeric value (the ICAT index), or
#'   \code{NA_real_} if computation is not possible.
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' if (max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE) > 0) {
#'   icat <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
#'   print(icat)
#' }
#'
#' @importFrom fastICA fastICA
#' @export
calc_icat <- function(patientID, tlsID, ldata = NULL) {

  ldata <- .resolve_ldata(ldata)
  .check_sample(patientID, ldata, required = c("x", "y", "tls_id_knn"))

  df <- ldata[[patientID]]

  tls_cells <- df[!is.na(df$tls_id_knn) & df$tls_id_knn == tlsID, ]

  if (nrow(tls_cells) < 3L) {
    message("calc_icat: TLS ", tlsID, " in '", patientID, "' has fewer than ",
            "3 cells. Returning NA.")
    return(NA_real_)
  }

  X  <- as.matrix(tls_cells[, c("x", "y")])
  mu <- colMeans(X, na.rm = TRUE)

  # Centre coordinates before ICA
  Xc        <- X
  Xc[, 1]   <- X[, 1] - mu[1]
  Xc[, 2]   <- X[, 2] - mu[2]

  ica_result <- tryCatch(
    fastICA::fastICA(Xc, n.comp = 2L, verbose = FALSE),
    error = function(e) {
      message("calc_icat: FastICA failed for TLS ", tlsID,
              " in '", patientID, "': ", conditionMessage(e),
              ". Returning NA.")
      NULL
    }
  )

  if (is.null(ica_result)) return(NA_real_)

  # Reconstruct data from ICA components and add back the mean
  Xhat <- ica_result$S %*% t(ica_result$A) + matrix(mu, nrow = nrow(X),
                                                     ncol = 2L, byrow = TRUE)

  # Marginal variances of reconstructed coordinates
  v <- diag(stats::cov(Xhat))    # length-2 vector: c(var_x, var_y)
  v1 <- v[1]; v2 <- v[2]

  # Trace-SD: always >= 0; equivalent to (sd_x + sd_y) scaled by sqrt
  trace_sd <- sqrt(v1 + v2 + 2 * sqrt(v1 * v2))

  # Normalise by number of cells so values are comparable across TLS sizes
  icat <- 100 * trace_sd / nrow(X)

  icat
}
