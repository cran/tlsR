#' Detect Tertiary Lymphoid Structures using a KNN-density approach
#'
#' Identifies TLS candidates by finding regions of high local B-cell density
#' that also contain a sufficient number of nearby T cells (B+T
#' co-localisation).  The algorithm proceeds in three steps:
#' \enumerate{
#'   \item For every B cell, compute the average inverse-kNN distance as a
#'         local density estimate.
#'   \item Retain B cells whose density exceeds \code{bcell_density_threshold}.
#'   \item Label connected components of dense B cells (DBSCAN-style, using
#'         \code{k}-NN graph edges) and keep components that also satisfy
#'         minimum B-cell count and T-cell proximity requirements.
#' }
#'
#' @param LSP Character. Sample name in \code{ldata}.
#' @param k Integer. Number of nearest neighbours used for density estimation
#'   (default \code{30}, calibrated for 0.325 um/px imaging).
#' @param bcell_density_threshold Numeric. Minimum average 1/k-distance (in
#'   microns) for a B cell to be considered locally dense (default \code{15}).
#' @param min_B_cells Integer. Minimum B cells per candidate TLS cluster
#'   (default \code{50}).
#' @param min_T_cells_nearby Integer. Minimum T cells within
#'   \code{max_distance_T} microns of the candidate cluster centre
#'   (default \code{30}).
#' @param max_distance_T Numeric. Search radius (microns) for T-cell proximity
#'   check (default \code{50}).
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#'@param expand_distance Numeric. The epanding radius from the boundary of the detected clusters to include other immune cells).
#' @return The input \code{ldata} list, with the data frame for \code{LSP}
#'   augmented by three new columns:
#'   \describe{
#'     \item{\code{tls_id_knn}}{Integer. \code{0} = non-TLS cell; positive
#'       integer = TLS cluster ID.}
#'     \item{\code{tls_center_x}}{Numeric. X coordinate of the TLS centre for
#'       TLS cells; \code{NA} otherwise.}
#'     \item{\code{tls_center_y}}{Numeric. Y coordinate of the TLS centre for
#'       TLS cells; \code{NA} otherwise.}
#'   }
#'
#' @examples
#' # Use a 70% sample of the data to keep CRAN check time under 10s.
#' # TLS detection requires sufficient cell density; 70% preserves
#' # the spatial structure needed for reliable detection.
#' # For production use, run on the full dataset (see \donttest{} below).
#' data(toy_ldata)
#' set.seed(42)
#' idx <- sample(nrow(toy_ldata[["ToySample"]]),
#'               size = floor(0.7 * nrow(toy_ldata[["ToySample"]])))
#' sub_ldata <- list(ToySample = toy_ldata[["ToySample"]][idx, ])
#' ldata <- detect_TLS("ToySample", k = 30, ldata = sub_ldata)
#' table(ldata[["ToySample"]]$tls_id_knn)
#' plot(ldata[["ToySample"]]$x, ldata[["ToySample"]]$y,
#'      col = ifelse(ldata[["ToySample"]]$tls_id_knn > 0, "red", "gray"),
#'      pch = 19, cex = 0.5, main = "Detected TLS (70% sample)")
#'
#' \donttest{
#' # Full dataset with default settings
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' table(ldata[["ToySample"]]$tls_id_knn)
#' plot(ldata[["ToySample"]]$x, ldata[["ToySample"]]$y,
#'      col = ifelse(ldata[["ToySample"]]$tls_id_knn > 0, "red", "gray"),
#'      pch = 19, cex = 0.5, main = "Detected TLS in toy data")
#' }
#'
#' @importFrom FNN get.knn
#' @export
detect_TLS <- function(LSP,
                       ldata,
                       k = 30L,
                       bcell_density_threshold = 10,
                       min_B_cells = 50L,
                       min_T_cells_nearby = 10L,
                       max_distance_T = 50,
                       expand_distance = 80) {
  
  # ----------------------------
  # Basic checks
  # ----------------------------
  if (is.null(ldata)) stop("ldata is NULL")
  if (!LSP %in% names(ldata)) stop("LSP not found in ldata")
  
  df <- ldata[[LSP]]
  
  if (!all(c("x", "y", "phenotype") %in% colnames(df))) {
    stop("Data must contain x, y, phenotype columns")
  }
  
  # ----------------------------
  # Initialize
  # ----------------------------
  df$tls_id_knn <- 0L
  df$tls_center_x <- NA_real_
  df$tls_center_y <- NA_real_
  
  # ----------------------------
  # Define B and T cells (SELF-CONTAINED)
  # ----------------------------
  b_idx <- which(grepl("B", df$phenotype, ignore.case = TRUE))
  t_idx <- which(grepl("T", df$phenotype, ignore.case = TRUE))
  
  if (length(b_idx) < k + 1L) {
    message("Not enough B cells for clustering")
    ldata[[LSP]] <- df
    return(ldata)
  }
  
  # ----------------------------
  # Coordinates
  # ----------------------------
  b_coords <- as.matrix(df[b_idx, c("x", "y")])
  t_coords <- if (length(t_idx) > 0)
    as.matrix(df[t_idx, c("x", "y")])
  else NULL
  
  # ----------------------------
  # Density estimate via kNN
  # ----------------------------
  knn_res <- FNN::get.knn(b_coords, k = k)
  avg_dist <- rowMeans(knn_res$nn.dist)
  density <- 1 / (avg_dist + .Machine$double.eps)
  
  dense_mask <- density >= (1 / bcell_density_threshold)
  
  if (sum(dense_mask) == 0) {
    message("No dense B-cell regions found")
    ldata[[LSP]] <- df
    return(ldata)
  }
  
  dense_b_coords <- b_coords[dense_mask, , drop = FALSE]
  dense_idx <- b_idx[dense_mask]
  
  # ----------------------------
  # Local kNN graph
  # ----------------------------
  k_local <- min(k, nrow(dense_b_coords) - 1L)
  if (k_local < 1L) {
    ldata[[LSP]] <- df
    return(ldata)
  }
  
  knn_dense <- FNN::get.knn(dense_b_coords, k = k_local)
  
  # ----------------------------
  # Connected components (CORE TLS)
  # ----------------------------
  component_id <- integer(nrow(dense_b_coords))
  tls_counter <- 0L
  
  for (i in seq_len(nrow(dense_b_coords))) {
    
    if (component_id[i] != 0L) next
    tls_counter <- tls_counter + 1L
    
    queue <- i
    
    while (length(queue) > 0) {
      cur <- queue[1]
      queue <- queue[-1]
      
      if (component_id[cur] != 0L) next
      
      component_id[cur] <- tls_counter
      
      nbrs <- knn_dense$nn.index[cur, ]
      nbrs <- nbrs[component_id[nbrs] == 0L]
      
      queue <- c(queue, nbrs)
    }
  }
  
  # ----------------------------
  # FILTER + ASSIGN CORE TLS
  # ----------------------------
  final_tls <- 0L
  
  for (cid in seq_len(tls_counter)) {
    
    members <- which(component_id == cid)
    
    if (length(members) < min_B_cells) next
    
    cx <- mean(dense_b_coords[members, 1])
    cy <- mean(dense_b_coords[members, 2])
    
    # T-cell validation
    n_t_near <- 0L
    if (!is.null(t_coords)) {
      dists_t <- sqrt((t_coords[,1] - cx)^2 + (t_coords[,2] - cy)^2)
      n_t_near <- sum(dists_t <= max_distance_T)
    }
    
    if (n_t_near < min_T_cells_nearby) next
    
    final_tls <- final_tls + 1L
    
    idx <- dense_idx[members]
    
    df$tls_id_knn[idx] <- final_tls
    df$tls_center_x[idx] <- cx
    df$tls_center_y[idx] <- cy
  }
  
  # ----------------------------
  # EXPANSION STEP (KEY ADDITION)
  # ----------------------------
  if (!is.null(expand_distance) && final_tls > 0) {
    
    for (cid in seq_len(final_tls)) {
      
      core_idx <- which(df$tls_id_knn == cid)
      if (length(core_idx) == 0) next
      
      cx <- mean(df$x[core_idx])
      cy <- mean(df$y[core_idx])
      
      dist_all <- sqrt((df$x - cx)^2 + (df$y - cy)^2)
      
      expand_idx <- which(dist_all <= expand_distance)
      
      df$tls_id_knn[expand_idx] <- cid
    }
  }
  
  message("Detected TLS: ", final_tls)
  
  ldata[[LSP]] <- df
  ldata
}