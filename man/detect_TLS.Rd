\name{detect_TLS}
\alias{detect_TLS}
\title{Detect Tertiary Lymphoid Structures using a KNN-density approach}
\usage{
detect_TLS(
  LSP,
  ldata,
  k = 30L,
  bcell_density_threshold = 10,
  min_B_cells = 50L,
  min_T_cells_nearby = 10L,
  max_distance_T = 50,
  expand_distance = 80
)
}
\arguments{
\item{LSP}{Character. Sample name in \code{ldata}.}

\item{k}{Integer. Number of nearest neighbours used for density estimation
(default \code{30}, calibrated for 0.325 um/px imaging).}

\item{bcell_density_threshold}{Numeric. Minimum average 1/k-distance (in
microns) for a B cell to be considered locally dense (default \code{15}).}

\item{min_B_cells}{Integer. Minimum B cells per candidate TLS cluster
(default \code{50}).}

\item{min_T_cells_nearby}{Integer. Minimum T cells within
\code{max_distance_T} microns of the candidate cluster centre
(default \code{30}).}

\item{max_distance_T}{Numeric. Search radius (microns) for T-cell proximity
check (default \code{50}).}

\item{expand_distance}{Integer. The extended values from the boundary of the deteced B-cells clusters that the T cells are bieng integrated (default \code{80}).}

\item{ldata}{Named list of data frames, or \code{NULL} to use the global
\code{ldata} object (deprecated; pass explicitly).}
}
\value{
The similarly formatted \code{ldata} list, with the data frame for \code{LSP}
augmented by three new columns:
\describe{
  \item{\code{tls_id_knn}}{Integer. \code{0} = non-TLS cell; positive
    integer = TLS cluster ID.}
  \item{\code{tls_center_x}}{Numeric. X coordinate of the TLS centre for
    TLS cells; \code{NA} otherwise.}
  \item{\code{tls_center_y}}{Numeric. Y coordinate of the TLS centre for
    TLS cells; \code{NA} otherwise.}
}
}
\description{
Identifies TLS candidates by finding regions of high local B-cell density
that also contain a sufficient number of nearby T cells (B+T
co-localisation). Phenotype labels \code{"B cell"} and \code{"B cells"}
(and their T-cell equivalents) are both accepted.
}
\examples{
# Use a 70\% sample of the data to keep CRAN check time under 10s.
# TLS detection requires sufficient cell density; 70\% preserves
# the spatial structure needed for reliable detection.
# For production use, run on the full dataset (see \donttest{} below).
data(toy_ldata)
set.seed(42)
idx <- sample(nrow(toy_ldata[["ToySample"]]),
              size = floor(0.7 * nrow(toy_ldata[["ToySample"]])))
sub_ldata <- list(ToySample = toy_ldata[["ToySample"]][idx, ])
ldata <- detect_TLS("ToySample", k = 30, ldata = sub_ldata)
table(ldata[["ToySample"]]$tls_id_knn)
plot(ldata[["ToySample"]]$x, ldata[["ToySample"]]$y,
     col = ifelse(ldata[["ToySample"]]$tls_id_knn > 0, "red", "gray"),
     pch = 19, cex = 0.5, main = "Detected TLS (70\% sample)")

\donttest{
# Full dataset with default settings
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
table(ldata[["ToySample"]]$tls_id_knn)
plot(ldata[["ToySample"]]$x, ldata[["ToySample"]]$y,
     col = ifelse(ldata[["ToySample"]]$tls_id_knn > 0, "red", "gray"),
     pch = 19, cex = 0.5, main = "Detected TLS in toy data")
}
}
