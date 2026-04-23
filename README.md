# tlsR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/tlsR)](https://CRAN.R-project.org/package=tlsR)
[![R-CMD-check](https://github.com/labsyspharm/tlsR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/labsyspharm/tlsR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**tlsR** provides fast, reproducible detection and spatial analysis of
Tertiary Lymphoid Structures (TLS) in multiplexed tissue imaging data
(mIHC, CODEX, IMC, etc.).

## Installation

```r
# Stable release from CRAN
install.packages("tlsR")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("labsyspharm/tlsR")
```

## Quick start

```r
library(tlsR)

data(toy_ldata)

# 1. Detect TLS
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)

# 2. Score each TLS with ICAT
calc_icat("ToySample", tlsID = 1, ldata = ldata)

# 3. Detect T-cell clusters outside TLS
ldata <- detect_tic("ToySample", ldata = ldata)

# 4. Tidy summary table
summarize_TLS(ldata)

# 5. Spatial plot
plot_TLS("ToySample", ldata = ldata)
```

## Data format

`ldata` is a **named list of data frames**, one per sample, each with columns:

| Column      | Type      | Description                     |
|-------------|-----------|---------------------------------|
| `x`         | numeric   | X coordinate (microns)          |
| `y`         | numeric   | Y coordinate (microns)          |
| `phenotype` | character | `"B cell"`, `"T cell"`, or other|

## Citation

Amiryousefi et al. (2025). *Detection and spatial analysis of tertiary
lymphoid structures in multiplexed tissue imaging.*
<https://doi.org/10.1101/2025.09.21.677465>

## License

MIT © Ali Amiryousefi
