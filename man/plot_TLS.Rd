\name{plot_TLS}
\alias{plot_TLS}
\title{Plot Spatial Map of TLS and T-cell Clusters}
\usage{
plot_TLS(
  sample,
  ldata = NULL,
  show_tic = TRUE,
  point_size = 0.5,
  alpha = 0.7,
  bg_alpha = 0.25,
  tic_size_mult = 1.8,
  tls_palette = c("#0072B2", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442"),
  tic_colour = "#E69F00",
  bg_colour = "grey80"
)
}
\arguments{
\item{sample}{Character. Sample name in \code{ldata}.}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}

\item{show_tic}{Logical. Colour T-cell clusters (if \code{detect_tic} has
been run) in a distinct colour?  Default \code{TRUE}.}

\item{point_size}{Numeric. Base point size for TLS cells and background
cells (default \code{0.5}).  TIC cells are drawn at
\code{point_size * tic_size_mult}.}

\item{alpha}{Numeric. Point transparency for TLS and TIC cells
(default \code{0.7}).}

\item{bg_alpha}{Numeric. Point transparency for background (non-TLS,
non-TIC) cells (default \code{0.25}).  Reducing this value pushes
background cells further behind the foreground structure.}

\item{tic_size_mult}{Numeric. Multiplier applied to \code{point_size} for
TIC cells so they appear slightly larger than background and TLS cells
(default \code{1.8}).}

\item{tls_palette}{Character vector of colours for TLS IDs.  Recycled if
there are more TLS than colours.  Default uses a colourblind-friendly
palette.}

\item{tic_colour}{Character. Colour for T-cell cluster cells
(default \code{"#E69F00"}).}

\item{bg_colour}{Character. Colour for non-TLS, non-TIC cells
(default \code{"grey80"}).}
}
\value{
A \code{ggplot} object (invisibly).  The plot is also printed unless
the return value is assigned.
}
\description{
Produces a \code{ggplot2} scatter plot of cell positions, coloured by
TLS membership, T-cell cluster membership, and background phenotype.

Background (non-TLS, non-TIC) cells are rendered with a lower alpha to
keep them visually recessive, while TIC cells are drawn slightly larger
than TLS cells so they stand out without dominating the plot.
}
\examples{
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
\donttest{
  p <- plot_TLS("ToySample", ldata = ldata)
}
}
