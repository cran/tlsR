\name{calc_icat}
\alias{calc_icat}
\title{Calculate ICAT (Immune Cell Arrangement Trace) Index}
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
A single numeric value (the ICAT index), or \code{NA_real_} if
computation is not possible (fewer than 3 cells, or FastICA did not
converge).
}
\description{
Quantifies the linear / organised spatial arrangement of cells within a
detected TLS by applying FastICA to the (x, y) coordinates of TLS cells
and returning the trace of the estimated mixing matrix. Higher values
indicate more structured (linear) cell organisation.

If the requested TLS contains fewer than 3 cells, or FastICA does not
converge, the function returns \code{NA_real_} with an informative message
rather than throwing an error.
}
\examples{
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
if (max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE) > 0) {
  icat <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
  icat
}
}
