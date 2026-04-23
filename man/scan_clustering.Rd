\name{scan_clustering}
\alias{scan_clustering}
\title{Scan Tissue for Local Immune Cell Clustering (Ripley's L Heatmap)}
\usage{
scan_clustering(
  ws = 500,
  sample,
  phenotype = c("T cells", "B cells", "Both"),
  plot = TRUE,
  creep = 1L,
  min_cells = 10L,
  min_phen_cells = 5L,
  label_cex = 1.1,
  ldata = NULL
)
}
\arguments{
\item{ws}{Numeric. Window side length in microns (default \code{500}).}

\item{sample}{Character. Sample name in \code{ldata}.}

\item{phenotype}{One of \code{"T cells"}, \code{"B cells"}, or
\code{"Both"}.}

\item{plot}{Logical. Draw the spatial clustering map?
(default \code{TRUE}).}

\item{creep}{Integer. Grid density factor; \code{creep = 2} overlaps
adjacent windows by half a window width, producing a smoother map
(default \code{1}).}

\item{min_cells}{Integer. Minimum total cell count required in a window
before it is analysed (default \code{10}).}

\item{min_phen_cells}{Integer. Minimum phenotype-specific cell count per
window (default \code{5}).}

\item{label_cex}{Numeric. Base character expansion for the CI numeric
labels drawn inside each window (default \code{1.1}).  Increase this
value if labels appear too small for your screen or output resolution.}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}
}
\value{
A named list with elements \code{B} and/or \code{T} (depending on
\code{phenotype}), each containing the \code{Lest} objects for all
qualifying windows of that phenotype.  Returned invisibly when
\code{plot = TRUE}.
}
\description{
Applies a sliding-window Ripley's L analysis across the tissue to produce
a spatial clustering map.  For each window a \emph{K-integral} index is
computed as the mean positive excess of the observed L function over its
theoretical CSR value.  When \code{plot = TRUE} a base-graphics spatial
map is drawn with LOESS-smoothed L-excess curves and numeric CI labels
overlaid inside each qualifying window, plus a legend identifying point
and curve colours.
}
\details{
The K-integral clustering index for window \eqn{w} is:
\deqn{\text{CI}_w = \frac{1}{N_+}\sum_{i:\,L_i > L_{\text{theo},i}}
      (L_i - L_{\text{theo},i})}
where \eqn{N_+} is the number of spatial lags where the observed L exceeds
the theoretical CSR value.

When \code{plot = TRUE} the map shows:
\itemize{
  \item All cells as small light-grey points.
  \item Phenotype cells (T cells green, B cells red).
  \item Navy dashed grid lines marking window boundaries.
  \item A LOESS-smoothed L-excess curve inside each qualifying window.
  \item A bold numeric CI label centred in the window.
  \item A legend identifying all point and curve colours.
}

When \code{phenotype = "Both"} two side-by-side panels are produced -
one for B cells and one for T cells - so the two clustering maps can be
compared directly on the same spatial layout.
}
\examples{
data(toy_ldata)
\donttest{
  L_models <- scan_clustering(
    ws        = 200,
    sample    = "ToySample",
    phenotype = "B cells",
    plot      = TRUE,
    ldata     = toy_ldata
  )
  cat("B-cell windows analysed:", length(L_models$B), "\n")

  # Side-by-side B and T cell panels
  L_both <- scan_clustering(
    ws        = 200,
    sample    = "ToySample",
    phenotype = "Both",
    plot      = TRUE,
    ldata     = toy_ldata
  )
}
}
