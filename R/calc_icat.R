#' Calculate ICAT (Immune Cell Arrangement Trace) Index
#'
#' Quantifies the linear / organised spatial arrangement of cells within a
#' detected TLS by applying FastICA to the (x, y) coordinates of TLS cells and
#' returning the trace of the estimated mixing matrix.  Higher values indicate
#' more structured (linear) cell organisation.
#'
#' @details
#' The ICAT index is defined as
#' \deqn{\text{ICAT} = \text{tr}(\mathbf{A})}
#' where \eqn{\mathbf{A}} is the 2x2 FastICA mixing matrix estimated from the
#' centred spatial coordinates of TLS cells.  A value near zero indicates
#' random arrangement; large positive values indicate high spatial linearity.
#'
#' If the requested TLS contains fewer than 3 cells, or FastICA does not
#' converge, the function returns \code{NA_real_} with an informative message
#' rather than throwing an error.
#'
#' @param patientID Character. Sample name in \code{ldata}.
#' @param tlsID Numeric or integer. TLS identifier (value of \code{tls_id_knn}
#'   for the TLS of interest).
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#'
#' @return A single numeric value (the ICAT index), or \code{NA_real_} if
#'   computation is not possible.
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' icat <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
#' icat
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

  coords <- as.matrix(tls_cells[, c("x", "y")])
  # Centre coordinates
  coords[, 1] <- coords[, 1] - mean(coords[, 1])
  coords[, 2] <- coords[, 2] - mean(coords[, 2])

  ica_result <- tryCatch(
    fastICA::fastICA(coords, n.comp = 2L, verbose = FALSE),
    error = function(e) {
      message("calc_icat: FastICA failed for TLS ", tlsID,
              " in '", patientID, "': ", conditionMessage(e),
              ". Returning NA.")
      NULL
    }
  )

  if (is.null(ica_result)) return(NA_real_)

  # Trace of the mixing matrix A
  A     <- ica_result$A
  icat  <- sum(diag(A))
  icat
}
