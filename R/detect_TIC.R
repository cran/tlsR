#' @name detect_tic
#' @title Detect Tumor-Infiltrating T-cell Clusters (TIC)
#' @description Uses HDBSCAN approximation on T cells outside TLS regions to find TIC.
#' @param sample Character. Sample name.
#' @param ldata Optional list.
#' @return Modified data frame for the sample (invisibly).
#' @export
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)  # Need TLS first
#' ldata <- detect_tic("ToySample", ldata = ldata)
#' table(ldata[["ToySample"]]$tcell_cluster_hdbscan)
utils::globalVariables(c("tls_id", "coarse_phen_vec", "row_index", "x", "y"))
detect_tic <- function(sample, ldata = NULL) {
  if (is.null(ldata)) {
    if (!exists("ldata", .GlobalEnv)) stop("ldata not found in global environment")
    ldata <- get("ldata", .GlobalEnv)
  }

  d <- ldata[[sample]]
  if (is.null(d)) stop("Sample not found in ldata")

  # Subset non-TLS T cells
  non_tls <- d$tls_id == 0
  tcells <- subset(d, coarse_phen_vec == "T cells" & non_tls)

  if (nrow(tcells) < 10) {
    message("Too few T cells outside TLS in ", sample)
    d$tcell_cluster_hdbscan <- -999
    d$tcell_dbscan <- 0
    ldata[[sample]] <- d
    return(invisible(d))
  }

  coords <- cbind(tcells$x, tcells$y)

  # Nearest neighbor stats
  nn <- RANN::nn2(coords, k = 2)$nn.dists[,2]
  message(sprintf("T cells (non-TLS) in %s: n=%d, median NN=%.2f um", sample, nrow(tcells), median(nn)))

  # DBSCAN as HDBSCAN approximation
  cl <- dbscan::dbscan(coords, eps = 50, minPts = 100)
  labels <- cl$cluster
  labels[labels == 0] <- -1  # noise

  d$tcell_cluster_hdbscan <- -999
  d$tcell_cluster_hdbscan[tcells$row_index] <- labels

  d$tcell_dbscan <- 0
  d$tcell_dbscan[non_tls] <- d$tcell_cluster_hdbscan[non_tls]

  ldata[[sample]] <- d
  return(invisible(d))
}
