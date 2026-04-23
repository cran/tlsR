#' Plot Spatial Map of TLS and T-cell Clusters
#'
#' Produces a \code{ggplot2} scatter plot of cell positions, coloured by
#' TLS membership, T-cell cluster membership, and background phenotype.
#'
#' Background (non-TLS, non-TIC) cells are rendered with a lower alpha to
#' keep them visually recessive.  TIC cells are drawn slightly larger than
#' background and TLS cells so they stand out without dominating the plot.
#'
#' @param sample Character. Sample name in \code{ldata}.
#' @param ldata Named list of data frames, or \code{NULL} to use the global
#'   \code{ldata} object (deprecated; pass explicitly).
#' @param show_tic Logical. Colour T-cell clusters (if \code{detect_tic} has
#'   been run) in a distinct colour?  Default \code{TRUE}.
#' @param point_size Numeric. Base point size passed to \code{geom_point} for
#'   TLS cells (default \code{0.5}).  Background cells use this same size;
#'   TIC cells are drawn at \code{point_size * tic_size_mult}.
#' @param alpha Numeric. Point transparency for TLS cells (default \code{0.7}).
#' @param bg_alpha Numeric. Point transparency for background cells
#'   (default \code{0.25}).  Lower values push background further behind
#'   the foreground structure.
#' @param tic_size_mult Numeric. Multiplier applied to \code{point_size} for
#'   TIC cells so they appear slightly larger (default \code{1.8}).
#' @param tls_palette Character vector of colours for TLS IDs.  Recycled if
#'   there are more TLS than colours.  Default uses a colourblind-friendly
#'   palette.
#' @param tic_colour Character. Colour for T-cell cluster cells
#'   (default \code{"#E69F00"}).
#' @param bg_colour Character. Colour for non-TLS, non-TIC cells
#'   (default \code{"grey80"}).
#'
#' @return A \code{ggplot} object (invisibly).  The plot is also printed unless
#'   the return value is assigned.
#'
#' @examples
#' data(toy_ldata)
#' ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#' \donttest{
#'   p <- plot_TLS("ToySample", ldata = ldata)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_classic labs
#'   scale_colour_manual theme element_text element_blank
#' @importFrom rlang .data
#' @export
plot_TLS <- function(sample,
                     ldata        = NULL,
                     show_tic     = TRUE,
                     point_size   = 0.5,
                     alpha        = 0.7,
                     bg_alpha     = 0.25,
                     tic_size_mult = 1.8,
                     tls_palette  = c("#0072B2", "#009E73", "#CC79A7",
                                      "#D55E00", "#56B4E9", "#F0E442"),
                     tic_colour   = "#E69F00",
                     bg_colour    = "grey80") {

  ldata <- .resolve_ldata(ldata)
  .check_sample(sample, ldata, required = c("x", "y", "tls_id_knn"))

  df <- ldata[[sample]]

  # -- Build grouping factor ------------------------------------------------
  tls_ids <- sort(unique(df$tls_id_knn[!is.na(df$tls_id_knn) &
                                         df$tls_id_knn > 0]))
  n_tls   <- length(tls_ids)
  tls_cols <- if (n_tls > 0) rep_len(tls_palette, n_tls) else character(0)

  has_tic <- show_tic && "tcell_cluster_hdbscan" %in% names(df)

  df$.group <- "Background"

  if (n_tls > 0) {
    for (i in seq_along(tls_ids)) {
      mask <- !is.na(df$tls_id_knn) & df$tls_id_knn == tls_ids[i]
      df$.group[mask] <- paste0("TLS ", tls_ids[i])
    }
  }

  if (has_tic) {
    tic_mask <- !is.na(df$tcell_cluster_hdbscan) &
                df$tcell_cluster_hdbscan > 0 &
                (is.na(df$tls_id_knn) | df$tls_id_knn == 0)
    df$.group[tic_mask] <- "TIC"
  }

  group_levels <- c("Background",
                    if (n_tls > 0) paste0("TLS ", tls_ids),
                    if (has_tic)   "TIC")
  df$.group <- factor(df$.group, levels = group_levels)

  colour_map <- c(Background = bg_colour)
  if (n_tls > 0)
    colour_map <- c(colour_map,
                    stats::setNames(tls_cols, paste0("TLS ", tls_ids)))
  if (has_tic)
    colour_map <- c(colour_map, TIC = tic_colour)

  # -- Per-group aesthetics (size and alpha) --------------------------------
  size_map  <- c(Background = point_size)
  alpha_map <- c(Background = bg_alpha)

  if (n_tls > 0) {
    tls_names <- paste0("TLS ", tls_ids)
    size_map[tls_names]  <- point_size
    alpha_map[tls_names] <- alpha
  }
  if (has_tic) {
    size_map["TIC"]  <- point_size * tic_size_mult
    alpha_map["TIC"] <- alpha
  }

  # Map each row to its size and alpha value
  df$.size  <- size_map[as.character(df$.group)]
  df$.alpha <- alpha_map[as.character(df$.group)]

  n_tic_cells <- if (has_tic) sum(df$.group == "TIC", na.rm = TRUE) else 0L
  subtitle_txt <- paste0(n_tls, " TLS detected",
                         if (has_tic && n_tic_cells > 0)
                           paste0("; ", n_tic_cells, " TIC cells") else "")

  # -- Build plot in two layers so size and alpha differ by group -----------
  # Layer 1: background (plotted first so it sits behind everything)
  df_bg  <- df[df$.group == "Background", ]
  df_fg  <- df[df$.group != "Background", ]

  p <- ggplot2::ggplot() +
    # Background - low alpha, small points
    ggplot2::geom_point(
      data  = df_bg,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.group),
      size  = point_size,
      alpha = bg_alpha
    ) +
    # Foreground TLS + TIC - higher alpha; TIC larger via mapped size
    ggplot2::geom_point(
      data  = df_fg,
      ggplot2::aes(x = .data$x, y = .data$y,
                   colour = .data$.group, size = .data$.size),
      alpha = alpha
    ) +
    ggplot2::scale_colour_manual(values = colour_map, name = "") +
    ggplot2::scale_size_identity() +   # use .size column directly
    ggplot2::theme_classic() +
    ggplot2::labs(
      title    = sample,
      subtitle = subtitle_txt,
      x        = "x (um)",
      y        = "y (um)"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40"),
      legend.key    = ggplot2::element_blank()
    )

  invisible(p)
}
