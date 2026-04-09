\name{plot_TLS}
\alias{plot_TLS}
\title{Plot Spatial Map of TLS and T-cell Clusters}
\usage{
plot_TLS(
  sample,
  ldata = NULL,
  show_tic = TRUE,
  point_size = 0.4,
  alpha = 0.6,
  tls_palette = c("#0072B2", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442"),
  tic_colour = "#E69F00",
  bg_colour = "grey80"
)
}
\arguments{
\item{sample}{Character. Sample name in \code{ldata}.}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}

\item{show_tic}{Logical. Colour T-cell clusters in a distinct colour?
Default \code{TRUE}.}

\item{point_size}{Numeric. Point size (default \code{0.4}).}

\item{alpha}{Numeric. Point transparency (default \code{0.6}).}

\item{tls_palette}{Character vector of colours for TLS IDs. Recycled if
there are more TLS than colours.}

\item{tic_colour}{Character. Colour for T-cell cluster cells
(default \code{"#E69F00"}).}

\item{bg_colour}{Character. Colour for background cells
(default \code{"grey80"}).}
}
\value{
A \code{ggplot} object (invisibly). Build and print the returned
object to display the plot.
}
\description{
Produces a \code{ggplot2} scatter plot of cell positions, coloured by
TLS membership, T-cell cluster membership, and background phenotype.
}
\examples{
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
\donttest{
  p <- plot_TLS("ToySample", ldata = ldata)
  print(p)
}
}
