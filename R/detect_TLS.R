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
#'
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
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' table(ldata[["ToySample"]]$tls_id_knn)
#' plot(ldata[["ToySample"]]$x, ldata[["ToySample"]]$y,
#'      col = ifelse(ldata[["ToySample"]]$tls_id_knn > 0, "red", "gray"),
#'      pch = 19, cex = 0.5, main = "Detected TLS in toy data")
#'
#' @importFrom FNN get.knn
#' @export
detect_TLS <- function(LSP,
                       k                       = 30L,
                       bcell_density_threshold = 10,
                       min_B_cells             = 50L,
                       min_T_cells_nearby      = 20L,
                       max_distance_T          = 50,
                       ldata                   = NULL) {

  ldata <- .resolve_ldata(ldata)
  .check_sample(LSP, ldata)

  df <- ldata[[LSP]]

  # Initialise output columns
  df$tls_id_knn   <- 0L
  df$tls_center_x <- NA_real_
  df$tls_center_y <- NA_real_

  # Subset B cells and T cells (accepts both "B cell" and "B cells")
  b_idx <- which(.is_bcell(df$phenotype))
  t_idx <- which(.is_tcell(df$phenotype))

  if (length(b_idx) < k + 1L) {
    message("Sample '", LSP, "': fewer B cells (", length(b_idx),
            ") than k+1 (", k + 1L, "); no TLS detected.")
    ldata[[LSP]] <- df
    return(ldata)
  }

  b_coords <- as.matrix(df[b_idx, c("x", "y")])
  t_coords <- if (length(t_idx) > 0) as.matrix(df[t_idx, c("x", "y")]) else NULL

  # Step 1: KNN density for B cells
  knn_res  <- FNN::get.knn(b_coords, k = k)
  avg_dist <- rowMeans(knn_res$nn.dist)
  density  <- 1 / (avg_dist + .Machine$double.eps)

  dense_mask <- density >= (1 / bcell_density_threshold)

  if (sum(dense_mask) == 0) {
    message("Sample '", LSP, "': no B cells exceed density threshold; ",
            "no TLS detected.")
    ldata[[LSP]] <- df
    return(ldata)
  }

  # Step 2: Connected-component labelling on dense B cells
  dense_b_coords <- b_coords[dense_mask, , drop = FALSE]
  n_dense        <- nrow(dense_b_coords)
  k_local        <- min(k, n_dense - 1L)

  if (k_local < 1L) {
    ldata[[LSP]] <- df
    return(ldata)
  }

  knn_dense    <- FNN::get.knn(dense_b_coords, k = k_local)
  component_id <- integer(n_dense)
  current_id   <- 0L

  for (i in seq_len(n_dense)) {
    if (component_id[i] != 0L) next
    current_id <- current_id + 1L
    queue <- i
    while (length(queue) > 0) {
      cur   <- queue[1]
      queue <- queue[-1]
      if (component_id[cur] != 0L) next
      component_id[cur] <- current_id
      nbrs  <- knn_dense$nn.index[cur, ]
      nbrs  <- nbrs[component_id[nbrs] == 0L]
      queue <- c(queue, nbrs)
    }
  }

  # Step 3: Apply size and T-cell proximity filters
  tls_counter <- 0L

  for (cid in seq_len(current_id)) {
    members <- which(component_id == cid)
    if (length(members) < min_B_cells) next

    cx <- mean(dense_b_coords[members, 1])
    cy <- mean(dense_b_coords[members, 2])

    n_t_near <- 0L
    if (!is.null(t_coords) && nrow(t_coords) > 0) {
      dists_t  <- sqrt((t_coords[, 1] - cx)^2 + (t_coords[, 2] - cy)^2)
      n_t_near <- sum(dists_t <= max_distance_T, na.rm = TRUE)
    }

    if (n_t_near < min_T_cells_nearby) next

    tls_counter <- tls_counter + 1L
    dense_b_original_idx <- b_idx[dense_mask][members]
    df$tls_id_knn[dense_b_original_idx]   <- tls_counter
    df$tls_center_x[dense_b_original_idx] <- cx
    df$tls_center_y[dense_b_original_idx] <- cy
  }

  if (tls_counter == 0L)
    message("Sample '", LSP, "': B-cell clusters found but none met T-cell ",
            "proximity threshold. Consider lowering `min_T_cells_nearby` or ",
            "increasing `max_distance_T`.")
  else
    message("Sample '", LSP, "': ", tls_counter, " TLS detected.")

  ldata[[LSP]] <- df
  ldata
}
