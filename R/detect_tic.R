#' Detect Tumor-Infiltrating T-cell Clusters (TIC)
#'
#' Applies HDBSCAN to T cells that lie \emph{outside} of previously detected
#' TLS regions to identify spatially compact T-cell clusters (TIC).
#'
#' @details
#' Only T cells with \code{tls_id_knn == 0} (i.e. not already assigned to a
#' TLS) are considered.  HDBSCAN (\code{\link[dbscan]{hdbscan}}) is run on
#' their (x, y) coordinates.  Cells assigned cluster label \code{0} by HDBSCAN
#' are treated as noise and receive \code{tcell_cluster_hdbscan = 0}.
#'
#' @param sample Character. Sample name in \code{ldata}.
#' @param min_pts Integer. HDBSCAN \code{minPts} parameter: minimum cluster
#'   size (default \code{10}).  Smaller values detect more, smaller clusters.
#' @param min_cluster_size Integer. Minimum number of T cells for a HDBSCAN
#'   cluster to be retained; smaller clusters are merged back into noise
#'   (label \code{0}).  Default \code{10}.
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#'
#' @return The input \code{ldata} list with the sample data frame augmented by
#'   one new column:
#'   \describe{
#'     \item{\code{tcell_cluster_hdbscan}}{Integer. \code{0} = noise / not a
#'       T-cell cluster; positive integer = TIC cluster ID.}
#'   }
#'   All non-T-cell rows receive \code{NA}.
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' ldata <- detect_tic("ToySample", ldata = ldata)
#' table(ldata[["ToySample"]]$tcell_cluster_hdbscan, useNA = "ifany")
#'
#' @importFrom dbscan hdbscan
#' @export
detect_tic <- function(sample,
                       min_pts          = 10L,
                       min_cluster_size = 10L,
                       ldata            = NULL) {

  ldata <- .resolve_ldata(ldata)
  .check_sample(sample, ldata, required = c("x", "y", "phenotype", "tls_id_knn"))

  df <- ldata[[sample]]

  df$tcell_cluster_hdbscan <- NA_integer_

  # T cells outside any TLS (accepts both "T cell" and "T cells")
  t_mask <- .is_tcell(df$phenotype) &
            !is.na(df$tls_id_knn) & df$tls_id_knn == 0L

  t_idx <- which(t_mask)

  if (length(t_idx) < min_pts) {
    message("detect_tic: sample '", sample, "' has fewer than `min_pts` (",
            min_pts, ") non-TLS T cells. Skipping HDBSCAN.")
    df$tcell_cluster_hdbscan[t_idx] <- 0L
    ldata[[sample]] <- df
    return(ldata)
  }

  t_coords <- as.matrix(df[t_idx, c("x", "y")])

  hdb <- tryCatch(
    dbscan::hdbscan(t_coords, minPts = min_pts),
    error = function(e) {
      message("detect_tic: HDBSCAN failed for '", sample, "': ",
              conditionMessage(e))
      NULL
    }
  )

  if (is.null(hdb)) {
    df$tcell_cluster_hdbscan[t_idx] <- 0L
    ldata[[sample]] <- df
    return(ldata)
  }

  labels <- hdb$cluster

  if (min_cluster_size > 1L) {
    cluster_sizes  <- table(labels[labels > 0L])
    small_clusters <- as.integer(
      names(cluster_sizes[cluster_sizes < min_cluster_size]))
    labels[labels %in% small_clusters] <- 0L
  }

  unique_ids <- sort(unique(labels[labels > 0L]))
  if (length(unique_ids) > 0) {
    mapping <- stats::setNames(seq_along(unique_ids), unique_ids)
    labels[labels > 0L] <- mapping[as.character(labels[labels > 0L])]
  }

  df$tcell_cluster_hdbscan[t_idx] <- as.integer(labels)

  n_tic <- length(unique(labels[labels > 0L]))
  if (n_tic > 0)
    message("detect_tic: ", n_tic, " T-cell cluster(s) detected in '",
            sample, "'.")
  else
    message("detect_tic: no significant T-cell clusters found in '",
            sample, "'.")

  ldata[[sample]] <- df
  ldata
}
