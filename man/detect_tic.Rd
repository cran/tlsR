\name{detect_tic}
\alias{detect_tic}
\title{Detect T-cell Immune Clusters (TIC) on the Tissue}
\usage{
detect_tic(sample, min_pts = 10L, min_cluster_size = 10L, ldata = NULL)
}
\arguments{
\item{sample}{Character. Sample name in \code{ldata}.}

\item{min_pts}{Integer. HDBSCAN \code{minPts} parameter: minimum cluster
size (default \code{10}). Smaller values detect more, smaller clusters.}

\item{min_cluster_size}{Integer. Minimum number of T cells for a HDBSCAN
cluster to be retained; smaller clusters are merged back into noise
(label \code{0}). Default \code{10}.}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}
}
\value{
The input \code{ldata} list with the sample data frame augmented by
one new column:
\describe{
  \item{\code{tcell_cluster_hdbscan}}{Integer. \code{0} = noise / not a
    T-cell cluster; positive integer = TIC cluster ID. Non-T-cell rows
    receive \code{NA}.}
}
}
\description{
Applies HDBSCAN to T cells that lie outside of previously detected TLS
regions to identify spatially compact T-cell clusters (TIC). Phenotype
labels \code{"T cell"} and \code{"T cells"} are both accepted.
}
\examples{
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
ldata <- detect_tic("ToySample", ldata = ldata)
table(ldata[["ToySample"]]$tcell_cluster_hdbscan, useNA = "ifany")
}
