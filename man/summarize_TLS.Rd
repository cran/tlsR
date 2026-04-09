\name{summarize_TLS}
\alias{summarize_TLS}
\title{Summarize Detected TLS Across Samples}
\usage{
summarize_TLS(ldata, calc_icat_scores = FALSE)
}
\arguments{
\item{ldata}{Named list of data frames as returned by \code{\link{detect_TLS}}
(and optionally \code{\link{detect_tic}}).}

\item{calc_icat_scores}{Logical. Should ICAT scores be computed for each TLS
and appended as a list-column? Default \code{FALSE}.}
}
\value{
A \code{data.frame} with columns:
\describe{
  \item{\code{sample}}{Sample name.}
  \item{\code{n_TLS}}{Number of TLS detected.}
  \item{\code{total_cells}}{Total cells in the sample.}
  \item{\code{TLS_cells}}{Number of cells assigned to any TLS.}
  \item{\code{TLS_fraction}}{Fraction of all cells that are TLS cells.}
  \item{\code{mean_TLS_size}}{Mean cells per TLS (\code{NA} if n_TLS = 0).}
  \item{\code{n_TIC}}{Number of T-cell clusters detected by
    \code{\link{detect_tic}} (\code{NA} if not yet run).}
  \item{\code{icat_scores}}{List-column of ICAT scores per TLS (only when
    \code{calc_icat_scores = TRUE}).}
}
}
\description{
Produces a tidy \code{data.frame} with one row per sample summarising the
number of detected TLS, their sizes, and (optionally) ICAT scores.
}
\examples{
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
summarize_TLS(ldata)
}
