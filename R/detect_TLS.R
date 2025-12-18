#' @name detect_TLS
#' @title Detect Tertiary Lymphoid Structures using a KNN-density approach
#' @description This function identifies TLS candidates germinated with B cells (BIC) by: 1. Finding regions of high local B-cell density (KNN-based). 2. Requiring a minimum number of neighbouring T cells (T-B cell co-localisation). 3. Applying sensible size and shape filters.
#' @param LSP Character. Sample name in the global \code{ldata} list.
#' @param k Integer. Number of nearest neighbours to consider (default 30 - works very well on 0.325 um/px imaging).
#' @param bcell_density_threshold Numeric. Minimum average 1/k-distance for B cells to be considered "dense" (default 15 um).
#' @param min_B_cells Integer. Minimum number of B cells in a candidate TLS (default 50).
#' @param min_T_cells_nearby Integer. Minimum T cells within 50 um of the B-cell cluster centre (default 30).
#' @param max_distance_T Numeric. Radius (um) to search for surrounding T cells (default 50).
#' @param ldata Optional. Named list of data frames. If \code{NULL}, uses global \code{ldata}.
#'
#' @return The original data frame with two new columns:
#' \item{tls_id_knn}{0 = non-TLS, positive integer = TLS cluster ID}
#' \item{tls_center_x, tls_center_y}{Coordinates of detected TLS centres (only for TLS cells)}
#'
#' @importFrom RANN nn2
#' @importFrom stats kmeans quantile
#' @export
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' table(ldata[["ToySample"]]$tls_id_knn)
#' plot(ldata[["ToySample"]]$x, ldata[["ToySample"]]$y,
#'      col = ifelse(ldata[["ToySample"]]$tls_id_knn > 0, "red", "gray"),
#'      pch = 19, cex = 0.5, main = "Detected TLS in toy data")
utils::globalVariables("coarse_phen_vec")
detect_TLS <- function(LSP,
                       k = 30,
                       bcell_density_threshold = 15,
                       min_B_cells = 50,
                       min_T_cells_nearby = 30,
                       max_distance_T = 50,
                       ldata = NULL) {

  if (is.null(ldata)) {
    if (!exists("ldata", envir = .GlobalEnv))
      stop("ldata not found in global environment.")
    ldata <- get("ldata", envir = .GlobalEnv)
  }

  d <- ldata[[LSP]]
  if (is.null(d)) stop(paste("Sample", LSP, "not found in ldata."))

  # Ensure clean starting point
  d$tls_id_knn <- 0
  d$tls_center_x <- NA_real_
  d$tls_center_y <- NA_real_

  # 1. Extract B cells only
  Bcells <- subset(d, coarse_phen_vec == "B cells")
  if (nrow(Bcells) == 0) {
    ldata[[LSP]] <- d
    message("No B cells found in ", LSP)
    return(invisible(ldata))
  }

  coords_B <- as.matrix(Bcells[, c("x", "y")])

  # 2. KNN distances for every B cell (k-th nearest B-cell neighbour)
  nn_B <- RANN::nn2(coords_B, k = k + 1)        # +1 because distance to itself = 0
  knn_dist_B <- nn_B$nn.dists[, k + 1]          # k-th nearest distance

  # 3. Local density score = 1 / mean distance to k nearest B cells to higher = denser
  density_score <- 1 / knn_dist_B

  # 4. Candidate B cells = very dense ones
  candidates <- which(knn_dist_B <= bcell_density_threshold & density_score > quantile(density_score, 0.85))

  if (length(candidates) == 0) {
    ldata[[LSP]] <- d
    message("No dense B-cell regions found in ", LSP)
    return(invisible(ldata))
  }

  candidate_coords <- coords_B[candidates, , drop = FALSE]

  # 5. Cluster the dense B-cell points (k-means is fast and excellent here)
  # Removed set.seed(42) as required by CRAN
  km <- kmeans(candidate_coords, centers = max(1, round(length(candidates)/50)), nstart = 10)
  cluster_labels <- km$cluster
  centers <- km$centers

  tls_id <- 1
  all_tls_cells <- integer(0)

  message(sprintf("Sample %s: Found %d candidate TLS cores (k=%d)", LSP, nrow(centers), k))

  # 6. Validate each candidate cluster
  for (i in seq_len(nrow(centers))) {
    centre_x <- centers[i, 1]
    centre_y <- centers[i, 2]

    # B cells belonging to this cluster
    cluster_B_idx_global <- Bcells$row_index[candidates[cluster_labels == i]]
    if (length(cluster_B_idx_global) < min_B_cells) next

    # T cells within max_distance_T um of this centre
    Tcells <- subset(d, coarse_phen_vec == "T cells")
    if (nrow(Tcells) == 0) next
    dist_to_centre <- sqrt((Tcells$x - centre_x)^2 + (Tcells$y - centre_y)^2)
    nearby_T <- sum(dist_to_centre <= max_distance_T)

    if (nearby_T >= min_T_cells_nearby) {
      # Accepted as real TLS
      d$tls_id_knn[cluster_B_idx_global] <- tls_id
      d$tls_center_x[cluster_B_idx_global] <- centre_x
      d$tls_center_y[cluster_B_idx_global] <- centre_y

      all_tls_cells <- c(all_tls_cells, cluster_B_idx_global)
      tls_id <- tls_id + 1
    }
  }

  # Optional: assign nearby T cells to the same TLS
  if (length(all_tls_cells) > 0) {
    tls_centres <- unique(d[all_tls_cells, c("tls_center_x", "tls_center_y")])
    Tcells <- subset(d, coarse_phen_vec == "T cells")
    for (j in seq_len(nrow(tls_centres))) {
      cx <- tls_centres$tls_center_x[j]
      cy <- tls_centres$tls_center_y[j]
      tid <- j
      close_T <- which(sqrt((Tcells$x - cx)^2 + (Tcells$y - cy)^2) <= max_distance_T)
      if (length(close_T) > 0) {
        d$tls_id_knn[Tcells$row_index[close_T]] <- tid
      }
    }
  }

  ldata[[LSP]] <- d
  message(sprintf("Detected %d TLS in %s", max(d$tls_id_knn), LSP))
  return(invisible(ldata))
}
