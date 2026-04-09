#' Scan Tissue for Local Immune Cell Clustering (Ripley's L)
#'
#' Applies a sliding-window centred L-function (CLF) analysis across the tissue
#' to identify spatially localised windows with statistically significant immune
#' cell clustering.  Each window is tested against a Monte Carlo envelope of
#' complete spatial randomness (CSR).
#'
#' @details
#' The tissue bounding box is divided into a regular grid of overlapping square
#' windows of side \code{ws} microns.  For each window that contains at least
#' \code{min_cells} cells of the requested phenotype, Ripley's L function is
#' estimated and an envelope is computed.  A window is declared "significant"
#' when the observed L curve exceeds the upper CSR envelope at any spatial lag.
#'
#' @param ws Numeric. Window side length in microns (default \code{500}).
#' @param sample Character. Sample name in \code{ldata}.
#' @param phenotype One of \code{"T cells"}, \code{"B cells"}, or \code{"Both"}.
#' @param plot Logical. Show a diagnostic plot for each significant window?
#'   (default \code{FALSE}).
#' @param creep Integer. Grid density factor; higher values give a finer
#'   scanning grid (default \code{1}).
#' @param nsim Integer. Number of Monte Carlo simulations for envelope
#'   estimation (default \code{39}, giving a pointwise significance level of
#'   0.05 under CSR).
#' @param min_cells Integer. Minimum cell count required in a window before it
#'   is analysed (default \code{10}).
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#'
#' @return A named list of results for significant windows.  Each element is
#'   itself a list with:
#'   \describe{
#'     \item{\code{Lest}}{The \code{spatstat} \code{Lest} object.}
#'     \item{\code{envelope}}{The Monte Carlo envelope object.}
#'     \item{\code{window_center}}{Numeric vector \code{c(cx, cy)} of window
#'       centre coordinates.}
#'     \item{\code{n_cells}}{Integer. Cell count in this window.}
#'   }
#'   Returns an empty list (invisibly) when no significant windows are found.
#'
#' @examples
#' data(toy_ldata)
#' \donttest{
#'   models <- scan_clustering(ws = 500, sample = "ToySample",
#'                             phenotype = "B cells", plot = FALSE,
#'                             nsim = 19, ldata = toy_ldata)
#'   length(models)
#' }
#'
#' @importFrom spatstat.geom ppp owin
#' @importFrom spatstat.explore Lest envelope
#' @export
scan_clustering <- function(ws        = 500,
                            sample,
                            phenotype = c("T cells", "B cells", "Both"),
                            plot      = FALSE,
                            creep     = 1L,
                            nsim      = 39L,
                            min_cells = 10L,
                            ldata     = NULL) {

  ldata     <- .resolve_ldata(ldata)
  phenotype <- match.arg(phenotype)
  .check_sample(sample, ldata)

  df <- ldata[[sample]]

  # Select cells (uses .is_bcell/.is_tcell to accept both singular and plural)
  sel <- switch(phenotype,
    "T cells" = .is_tcell(df$phenotype),
    "B cells" = .is_bcell(df$phenotype),
    "Both"    = .is_tcell(df$phenotype) | .is_bcell(df$phenotype)
  )

  sub <- df[sel, ]

  if (nrow(sub) == 0) {
    message("scan_clustering: no cells of phenotype '", phenotype,
            "' in sample '", sample, "'.")
    return(invisible(list()))
  }

  x_range <- range(sub$x, na.rm = TRUE)
  y_range <- range(sub$y, na.rm = TRUE)

  # Guard: tissue extent must be larger than window on both axes
  if ((x_range[2] - x_range[1]) <= ws || (y_range[2] - y_range[1]) <= ws) {
    message("scan_clustering: tissue extent smaller than window size; ",
            "try a smaller `ws`.")
    return(invisible(list()))
  }

  step   <- ws / creep
  x_ctrs <- seq(x_range[1] + ws / 2, x_range[2] - ws / 2, by = step)
  y_ctrs <- seq(y_range[1] + ws / 2, y_range[2] - ws / 2, by = step)

  results <- list()

  for (cx in x_ctrs) {
    for (cy in y_ctrs) {
      win_cells <- sub[
        sub$x >= cx - ws / 2 & sub$x <= cx + ws / 2 &
        sub$y >= cy - ws / 2 & sub$y <= cy + ws / 2, ]

      if (nrow(win_cells) < min_cells) next

      win <- tryCatch(
        spatstat.geom::owin(
          xrange = c(cx - ws / 2, cx + ws / 2),
          yrange = c(cy - ws / 2, cy + ws / 2)
        ),
        error = function(e) NULL
      )
      if (is.null(win)) next

      pp <- tryCatch(
        spatstat.geom::ppp(win_cells$x, win_cells$y, window = win),
        error = function(e) NULL
      )
      if (is.null(pp)) next

      L_obs <- tryCatch(
        spatstat.explore::Lest(pp, correction = "border"),
        error = function(e) NULL
      )
      if (is.null(L_obs)) next

      env <- tryCatch(
        spatstat.explore::envelope(pp, fun = spatstat.explore::Lest,
                                   nsim = nsim, verbose = FALSE,
                                   correction = "border"),
        error = function(e) NULL
      )
      if (is.null(env)) next

      is_sig <- any(env$obs > env$hi, na.rm = TRUE)

      if (is_sig) {
        key <- paste0("x", round(cx), "_y", round(cy))
        results[[key]] <- list(
          Lest          = L_obs,
          envelope      = env,
          window_center = c(cx, cy),
          n_cells       = nrow(win_cells)
        )
        if (plot) {
          graphics::plot(env,
            main = paste0(sample, " | centre (", round(cx), ", ", round(cy),
                          ") | n=", nrow(win_cells)))
        }
      }
    }
  }

  if (length(results) == 0)
    message("scan_clustering: no significant windows found in '", sample, "'.")
  else
    message("scan_clustering: ", length(results), " significant window(s) ",
            "found in '", sample, "'.")

  results
}
