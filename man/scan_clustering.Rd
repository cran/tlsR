\name{scan_clustering}
\alias{scan_clustering}
\title{Scan Tissue for Local Immune Cell Clustering (Ripley's L)}
\usage{
scan_clustering(
  ws = 500,
  sample,
  phenotype = c("T cells", "B cells", "Both"),
  plot = FALSE,
  creep = 1L,
  nsim = 39L,
  min_cells = 10L,
  ldata = NULL
)
}
\arguments{
\item{ws}{Numeric. Window side length in microns (default \code{500}).}

\item{sample}{Character. Sample name in \code{ldata}.}

\item{phenotype}{One of \code{"T cells"}, \code{"B cells"}, or \code{"Both"}.}

\item{plot}{Logical. Show a diagnostic plot for each significant window?
(default \code{FALSE}).}

\item{creep}{Integer. Grid density factor; higher values give a finer
scanning grid (default \code{1}).}

\item{nsim}{Integer. Number of Monte Carlo simulations for envelope
estimation (default \code{39}, giving a pointwise significance level of
0.05 under CSR).}

\item{min_cells}{Integer. Minimum cell count required in a window before it
is analysed (default \code{10}).}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}
}
\value{
A named list of results for significant windows. Each element is
itself a list with:
\describe{
  \item{\code{Lest}}{The \code{spatstat} \code{Lest} object.}
  \item{\code{envelope}}{The Monte Carlo envelope object.}
  \item{\code{window_center}}{Numeric vector \code{c(cx, cy)} of window
    centre coordinates.}
  \item{\code{n_cells}}{Integer. Cell count in this window.}
}
Returns an empty list (invisibly) when no significant windows are found.
}
\description{
Applies a sliding-window centred L-function analysis across the tissue to
identify spatially localised windows with statistically significant immune
cell clustering. Each window is tested against a Monte Carlo envelope of
complete spatial randomness (CSR).
}
\examples{
data(toy_ldata)
\donttest{
  models <- scan_clustering(ws = 500, sample = "ToySample",
                            phenotype = "B cells", plot = FALSE,
                            nsim = 19, ldata = toy_ldata)
  length(models)
}
}
