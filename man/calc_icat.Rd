\name{calc_icat}
\alias{calc_icat}
\title{Calculate ICAT (Independent Component Analysis Trace) Index for  TLS}
\usage{
calc_icat(patientID, tlsID, ldata = NULL)
}
\arguments{
\item{patientID}{Character. Sample name in \code{ldata}.}

\item{tlsID}{Numeric or integer. TLS identifier (value of \code{tls_id_knn}
for the TLS of interest).}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}
}
\value{
A single non-negative numeric value (the ICAT index), or
\code{NA_real_} if computation is not possible (fewer than 3 cells, or
FastICA did not converge).
}
\description{
Quantifies the spatial spread and linear organisation of cells within a
detected TLS.  FastICA is applied to the (x, y) coordinates of TLS cells
to estimate independent components; the mixing matrix is used to
reconstruct the data, and the ICAT index is defined as the normalised
trace-standard-deviation of the reconstructed coordinates.

The index is always non-negative because it measures the average
spatial spread per cell rather than the signed trace of the mixing
matrix (which can be negative due to ICA sign ambiguity).  Higher values
indicate a more spatially extended, structured cluster.
}
\details{
The ICAT index is computed as follows:
\enumerate{
  \item Centre the (x, y) coordinates of TLS cells.
  \item Run \code{fastICA} with 2 components.
  \item Reconstruct \eqn{\hat{X} = S A^T + \mu}.
  \item Let \eqn{v_1, v_2} be the marginal variances of \eqn{\hat{X}}.
  \item \deqn{\text{ICAT} = 100 \times
        \frac{\sqrt{v_1 + v_2 + 2\sqrt{v_1 v_2}}}{\text{nrow}(X)}}
}

If the requested TLS contains fewer than 3 cells, or FastICA does not
converge, the function returns \code{NA_real_} with an informative message
rather than throwing an error.
}
\examples{
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
if (max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE) > 0) {
  icat <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
  print(icat)
}
}
