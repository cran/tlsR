#' Scan Tissue for Local Immune Cell Clustering (Ripley's L Heatmap)
#'
#' Applies a sliding-window Ripley's L analysis across the tissue to map the
#' spatial intensity of immune cell clustering.  For each window the mean
#' positive excess of the observed L function over its theoretical CSR value
#' (the \emph{K-integral} CI) is computed and, when \code{plot = TRUE},
#' overlaid on a spatial scatter plot as LOESS-smoothed curves and numeric
#' labels.
#'
#' @details
#' The K-integral clustering index for window \eqn{w} is:
#' \deqn{\text{CI}_w = \frac{1}{N_+}\sum_{i:\,L_i > L_{\text{theo},i}}
#'       (L_i - L_{\text{theo},i})}
#' where \eqn{N_+} is the number of spatial lags where the observed L exceeds
#' the theoretical CSR value.
#'
#' When \code{plot = TRUE} the map shows:
#' \itemize{
#'   \item All cells as small light-grey points.
#'   \item Phenotype cells (T cells green, B cells red).
#'   \item Navy dashed grid lines marking window boundaries.
#'   \item A LOESS-smoothed L-excess curve inside each qualifying window.
#'   \item A numeric CI label centred in the window.
#'   \item A legend identifying point colours and curve colours.
#' }
#' When \code{phenotype = "Both"} two side-by-side panels are drawn  -  one for
#' B cells and one for T cells  -  so the clustering maps can be compared
#' directly.
#'
#' @param ws Numeric. Window side length in microns (default \code{500}).
#' @param sample Character. Sample name in \code{ldata}.
#' @param phenotype One of \code{"T cells"}, \code{"B cells"}, or
#'   \code{"Both"}.
#' @param plot Logical. Draw the spatial clustering map?  (default \code{TRUE}).
#' @param creep Integer. Grid density factor; \code{creep = 2} overlaps
#'   adjacent windows by half a window width, giving a smoother map
#'   (default \code{1}).
#' @param min_cells Integer. Minimum total cell count per window before it is
#'   analysed (default \code{10}).
#' @param min_phen_cells Integer. Minimum phenotype-specific cell count per
#'   window (default \code{5}).
#' @param label_cex Numeric. Base character expansion for the CI numeric labels
#'   drawn inside each window (default \code{1.1}).  Increase for larger text.
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#'
#' @return A named list with elements \code{B} and/or \code{T} (depending on
#'   \code{phenotype}), each containing the \code{Lest} objects for all
#'   qualifying windows of that phenotype.  Returned invisibly when
#'   \code{plot = TRUE}.
#'
#' @examples
#' data(toy_ldata)
#' \donttest{
#'   L_models <- scan_clustering(
#'     ws        = 200,
#'     sample    = "ToySample",
#'     phenotype = "B cells",
#'     plot      = TRUE,
#'     ldata     = toy_ldata
#'   )
#'   cat("B-cell windows analysed:", length(L_models$B), "\n")
#' }
#'
#' @importFrom spatstat.geom ppp owin
#' @importFrom spatstat.explore Lest
#' @importFrom graphics plot points abline lines text legend par mtext
#' @importFrom grDevices adjustcolor
#' @export
scan_clustering <- function(ws             = 500,
                            sample,
                            phenotype      = c("T cells", "B cells", "Both"),
                            plot           = TRUE,
                            creep          = 1L,
                            min_cells      = 10L,
                            min_phen_cells = 5L,
                            label_cex      = 1.1,
                            ldata          = NULL) {

  ldata     <- .resolve_ldata(ldata)
  phenotype <- match.arg(phenotype)
  .check_sample(sample, ldata)

  d   <- ldata[[sample]]
  pws <- ws / 0.65   # micron -> pixel window (assumes 0.65 um/px)

  # -- Internal: draw one spatial panel --------------------------------------
  .draw_panel <- function(d, pws, ws, sample, phen_label,
                          pt_col, line_col, text_col, line_lty,
                          panel_title) {

    graphics::plot(
      d$x, d$y,
      pch     = 19, cex = 0.01,
      cex.axis = 1, cex.lab = 1, cex.main = 1.3,
      col     = "lightgrey",
      main    = panel_title,
      xlab    = paste("Window size =", ws, "um"),
      ylab    = "",
      ylim    = range(d$y),
      xlim    = range(d$x),
      col.main = "navy"
    )

    # Tumour / epithelial cells (optional column)
    if ("panCKp" %in% names(d) && sum(d$panCKp, na.rm = TRUE) > 0) {
      ss <- d[d$panCKp == 1, ]
      graphics::points(ss$x, ss$y,
                       col = grDevices::adjustcolor("grey57", alpha.f = 0.2),
                       pch = 19, cex = 0.005)
    }

    # Phenotype cells
    phen_rows <- if (phen_label == "T cells")
                   d[.is_tcell(d$phenotype), ]
                 else
                   d[.is_bcell(d$phenotype), ]

    if (nrow(phen_rows) > 0)
      graphics::points(phen_rows$x, phen_rows$y,
                       col = grDevices::adjustcolor(pt_col, alpha.f = 0.4),
                       pch = 19, cex = 0.005)

    # Grid lines
    nx_g <- ceiling(max(d$x) / pws)
    ny_g <- ceiling(max(d$y) / pws)

    for (i in 0:ny_g)
      graphics::abline(h   = 1 + pws * i,
                       col = grDevices::adjustcolor("navy", alpha.f = 0.5),
                       lty = 2, lwd = 1.5)
    for (i in 0:ny_g)
      graphics::abline(h   = 1 + pws * i + pws / 2,
                       col = grDevices::adjustcolor("darkolivegreen", alpha.f = 0.5),
                       lty = 3, lwd = 1)
    for (j in 0:nx_g)
      graphics::abline(v   = 1 + pws * j,
                       col = grDevices::adjustcolor("navy", alpha.f = 0.5),
                       lty = 2, lwd = 1.5)

    # Legend for this panel
    legend_labels <- c("All cells", phen_label, "L-excess curve")
    legend_cols   <- c("lightgrey", pt_col, line_col)
    legend_lty    <- c(NA, NA, line_lty)
    legend_pch    <- c(19, 19, NA)
    legend_lwd    <- c(NA, NA, 2)

    if ("panCKp" %in% names(d) && sum(d$panCKp, na.rm = TRUE) > 0) {
      legend_labels <- c(legend_labels, "Tumour (panCK+)")
      legend_cols   <- c(legend_cols, "grey57")
      legend_lty    <- c(legend_lty, NA)
      legend_pch    <- c(legend_pch, 19)
      legend_lwd    <- c(legend_lwd, NA)
    }

    graphics::legend(
      "bottomleft",
      legend  = legend_labels,
      col     = legend_cols,
      lty     = legend_lty,
      pch     = legend_pch,
      lwd     = legend_lwd,
      pt.cex  = 1.4,
      cex     = 0.85,
      bty     = "n",
      bg      = grDevices::adjustcolor("white", alpha.f = 0.7)
    )
  }

  # -- Internal: compute L and overlay curves + labels on existing plot -------
  .run_phen <- function(d, phen_label, line_col, text_col,
                        text_y_frac, line_lty, pws, ws, creep,
                        min_cells, min_phen_cells, label_cex, do_plot) {

    nx <- ceiling(max(d$x) / pws)
    ny <- ceiling(max(d$y) / pws)

    L_models <- list()
    c_idx    <- 1L

    for (k in seq_len(creep)) {
      xstart <- k * round(pws / creep) - round(pws / creep)
      for (l in seq_len(creep)) {
        ystart <- l * round(pws / creep) - round(pws / creep)

        for (i in seq_len(nx)) {
          for (j in seq_len(ny)) {

            win_data <- d[
              d$x > xstart + pws * (i - 1) & d$x < xstart + pws * i &
              d$y > ystart + pws * (j - 1) & d$y < ystart + pws * j, ]

            if (nrow(win_data) < min_cells) next

            phen_data <- switch(phen_label,
              "T cells" = win_data[.is_tcell(win_data$phenotype), ],
              "B cells" = win_data[.is_bcell(win_data$phenotype), ]
            )

            if (nrow(phen_data) < min_phen_cells)              next
            if (diff(range(phen_data$x)) == 0 ||
                diff(range(phen_data$y)) == 0)                 next

            pp <- tryCatch(
              spatstat.geom::ppp(
                phen_data$x, phen_data$y,
                window = spatstat.geom::owin(
                  c(min(phen_data$x), max(phen_data$x)),
                  c(min(phen_data$y), max(phen_data$y))
                )
              ),
              error = function(e) NULL
            )
            if (is.null(pp)) next

            L <- tryCatch(
              spatstat.explore::Lest(pp, rmax = ws, correction = "border"),
              error = function(e) NULL
            )
            if (is.null(L)) next

            L_models[[c_idx]] <- L
            c_idx <- c_idx + 1L

            if (do_plot && k + l < 3L) {

              differences <- L$border - L$theo
              valid        <- is.finite(differences)

              # LOESS curve
              if (sum(valid) > 10L) {
                xvals <- seq_along(differences)
                sm    <- tryCatch(
                  stats::loess(differences[valid] ~ xvals[valid], span = 0.3),
                  error = function(e) NULL
                )
                if (!is.null(sm)) {
                  sm_vals <- stats::predict(sm, xvals)
                  sm_vals[!is.finite(sm_vals)] <- NA_real_
                  graphics::lines(
                    seq(xstart + pws * (i - 1), xstart + pws * i,
                        length.out = length(sm_vals)),
                    sm_vals + pws / 2 + pws * (j - 1),
                    col = line_col,
                    lwd = max(1, 2 * ws / 1500),
                    lty = line_lty
                  )
                }
              }

              # CI numeric label  -  use label_cex as the base, scaled mildly
              # by the CI value so high-clustering windows stand out slightly,
              # but never smaller than label_cex * 0.7.
              D  <- differences[is.finite(differences) & differences > 0]
              ci <- if (length(D) > 0L) mean(D, na.rm = TRUE) else 0
              scale_factor <- 1 + 0.3 * min(ci / max(1, ci + 500), 1)
              cex_use <- max(label_cex * 0.7, label_cex * scale_factor)

              graphics::text(
                x      = xstart + pws * (i - 0.5),
                y      = ystart + pws * text_y_frac + pws * (j - 1),
                labels = round(ci),
                col    = text_col,
                cex    = cex_use,
                font   = 2L   # bold so labels read clearly at any size
              )
            }
          }
        }
      }
    }

    L_models
  }

  # -- Layout: two panels for "Both", one panel otherwise --------------------
  old_par <- NULL
  if (plot && phenotype == "Both") {
    old_par <- graphics::par(mfrow = c(1, 2),
                             mar   = c(4, 3, 3, 1),
                             oma   = c(0, 0, 2, 0))
  }

  # -- B-cell pass -----------------------------------------------------------
  L_B <- list()
  if (phenotype %in% c("B cells", "Both")) {

    if (plot) {
      panel_title <- if (phenotype == "Both")
                       paste0(sample, "  -  B cells")
                     else
                       sample
      .draw_panel(d, pws, ws, sample,
                  phen_label  = "B cells",
                  pt_col      = "red",
                  line_col    = "plum1",
                  text_col    = "plum4",
                  line_lty    = 1L,
                  panel_title = panel_title)
    }

    L_B <- .run_phen(
      d              = d,
      phen_label     = "B cells",
      line_col       = "plum1",
      text_col       = "plum4",
      text_y_frac    = 0.5,
      line_lty       = 1L,
      pws            = pws,
      ws             = ws,
      creep          = creep,
      min_cells      = min_cells,
      min_phen_cells = min_phen_cells,
      label_cex      = label_cex,
      do_plot        = plot
    )

    message("scan_clustering [B cells]: ", length(L_B),
            " window(s) analysed in '", sample, "'.")
  }

  # -- T-cell pass -----------------------------------------------------------
  L_T <- list()
  if (phenotype %in% c("T cells", "Both")) {

    if (plot) {
      panel_title <- if (phenotype == "Both")
                       paste0(sample, "  -  T cells")
                     else
                       sample
      .draw_panel(d, pws, ws, sample,
                  phen_label  = "T cells",
                  pt_col      = "green",
                  line_col    = "darkgreen",
                  text_col    = "darkgreen",
                  line_lty    = 2L,
                  panel_title = panel_title)
    }

    L_T <- .run_phen(
      d              = d,
      phen_label     = "T cells",
      line_col       = "darkgreen",
      text_col       = "darkgreen",
      text_y_frac    = 0.5,
      line_lty       = 2L,
      pws            = pws,
      ws             = ws,
      creep          = creep,
      min_cells      = min_cells,
      min_phen_cells = min_phen_cells,
      label_cex      = label_cex,
      do_plot        = plot
    )

    message("scan_clustering [T cells]: ", length(L_T),
            " window(s) analysed in '", sample, "'.")
  }

  # -- "Both" super-title ----------------------------------------------------
  if (plot && phenotype == "Both") {
    graphics::mtext(paste("Clustering scan  - ", sample),
                    outer = TRUE, cex = 1.2, font = 2, col = "navy")
    graphics::par(old_par)   # restore layout
  }

  # -- Return ----------------------------------------------------------------
  result <- list()
  if (length(L_B) > 0L) result$B <- L_B
  if (length(L_T) > 0L) result$T <- L_T

  if (length(result) == 0L)
    message("scan_clustering: no windows met the cell-count thresholds in '",
            sample, "'.")

  invisible(result)
}
