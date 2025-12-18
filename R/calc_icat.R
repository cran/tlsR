#' @name calc_icat
#' @title Calculate ICAT (Immune Cell Arrangement Trace) Index
#' @description Quantifies linear/organized arrangement of cells within a TLS using FastICA.
#' @param patientID Character. Name of the sample in `ldata`.
#' @param tlsID Numeric/integer. TLS identifier.
#' @param ldata Named list of data frames (optional; defaults to global `ldata`).
#' @return Numeric ICAT value.
#' @export
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)  # First detect TLS
#' if (max(ldata[["ToySample"]]$tls_id_knn) > 0) {
#'   icat <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
#'   icat
#' }
utils::globalVariables(c("tls_id", "x", "y"))
calc_icat <- function(patientID, tlsID, ldata = NULL) {
  if (is.null(ldata)) {
    if (!exists("ldata", .GlobalEnv)) stop("ldata not found in global environment")
    ldata <- get("ldata", .GlobalEnv)
  }
  d <- ldata[[patientID]]
  if (is.null(d)) stop("Patient '", patientID, "' not found in ldata")

  df_xy <- subset(d, tls_id == tlsID)
  if (nrow(df_xy) == 0) stop("No cells for this tlsID")

  X <- as.matrix(df_xy[, c("x", "y")] * 2)
  mu <- colMeans(X, na.rm = TRUE)

  ica <- fastICA::fastICA(X, n.comp = 2)
  Xhat <- ica$S %*% t(ica$A) + mu

  traceSD <- sqrt(diag(cov(Xhat))[1] + diag(cov(Xhat))[2] +
                    2 * sqrt(diag(cov(Xhat))[1] * diag(cov(Xhat))[2]))
  ICAT <- 100 * traceSD / nrow(X)
  return(ICAT)
}
