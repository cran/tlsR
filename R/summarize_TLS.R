#' Summarize Detected TLS Across Samples
#'
#' Produces a tidy \code{data.frame} with one row per sample summarising the
#' number of detected TLS, their sizes, and (optionally) ICAT scores.
#'
#' @param ldata Named list of data frames as returned by \code{\link{detect_TLS}}
#'   (and optionally \code{\link{detect_tic}}).
#' @param calc_icat_scores Logical. Should ICAT scores be computed for each TLS
#'   and appended as a list-column?  Set to \code{FALSE} to skip FastICA
#'   (faster for large datasets).  Default \code{FALSE}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{sample}}{Sample name.}
#'     \item{\code{n_TLS}}{Number of TLS detected.}
#'     \item{\code{total_cells}}{Total cells in the sample.}
#'     \item{\code{TLS_cells}}{Number of cells assigned to any TLS.}
#'     \item{\code{TLS_fraction}}{Fraction of all cells that are TLS cells.}
#'     \item{\code{mean_TLS_size}}{Mean cells per TLS (\code{NA} if n_TLS = 0).}
#'     \item{\code{n_TIC}}{Number of T-cell clusters (if \code{detect_tic} has
#'       been run; \code{NA} otherwise).}
#'     \item{\code{icat_scores}}{List-column of numeric ICAT scores per TLS
#'       (only when \code{calc_icat_scores = TRUE}).}
#'   }
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' summarize_TLS(ldata)
#'
#' @export
summarize_TLS <- function(ldata, calc_icat_scores = FALSE) {

  ldata <- .resolve_ldata(ldata)

  if (length(ldata) == 0)
    stop("`ldata` is empty.", call. = FALSE)

  rows <- lapply(names(ldata), function(s) {
    df <- ldata[[s]]

    if (!"tls_id_knn" %in% names(df)) {
      warning("Sample '", s, "' does not have a `tls_id_knn` column. ",
              "Run detect_TLS() first.", call. = FALSE)
      return(data.frame(
        sample       = s,
        n_TLS        = NA_integer_,
        total_cells  = nrow(df),
        TLS_cells    = NA_integer_,
        TLS_fraction = NA_real_,
        mean_TLS_size = NA_real_,
        n_TIC        = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    tls_ids       <- sort(unique(df$tls_id_knn[df$tls_id_knn > 0]))
    n_tls         <- length(tls_ids)
    tls_cell_cnt  <- sum(df$tls_id_knn > 0, na.rm = TRUE)
    tls_sizes     <- if (n_tls > 0)
                       vapply(tls_ids, function(id) sum(df$tls_id_knn == id), integer(1L))
                     else integer(0)

    n_tic <- if ("tcell_cluster_hdbscan" %in% names(df)) {
      length(unique(df$tcell_cluster_hdbscan[
        !is.na(df$tcell_cluster_hdbscan) & df$tcell_cluster_hdbscan > 0]))
    } else NA_integer_

    result <- data.frame(
      sample        = s,
      n_TLS         = n_tls,
      total_cells   = nrow(df),
      TLS_cells     = tls_cell_cnt,
      TLS_fraction  = tls_cell_cnt / nrow(df),
      mean_TLS_size = if (n_tls > 0) mean(tls_sizes) else NA_real_,
      n_TIC         = n_tic,
      stringsAsFactors = FALSE
    )

    if (calc_icat_scores && n_tls > 0) {
      scores <- vapply(tls_ids,
        function(id) calc_icat(s, tlsID = id, ldata = ldata),
        numeric(1L))
      result$icat_scores <- list(stats::setNames(scores, paste0("TLS", tls_ids)))
    }

    result
  })

  do.call(rbind, rows)
}
