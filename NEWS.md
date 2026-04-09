# tlsR 0.2.0

## New functions

* `summarize_TLS()` — produces a tidy `data.frame` with one row per sample
  summarising detected TLS count, size, cell fractions, T-cell clusters, and
  (optionally) per-TLS ICAT scores.

* `plot_TLS()` — ggplot2-based spatial scatter plot of TLS and T-cell cluster
  membership, using a colourblind-friendly palette. Returns the `ggplot`
  object invisibly for further customisation.

## Improvements to existing functions

* `calc_icat()` — now returns `NA_real_` with an informative `message()` when
  a TLS has fewer than 3 cells or when FastICA fails to converge, instead of
  throwing an error.

* `detect_tic()` — now returns the updated `ldata` list visibly (was
  invisible). Two new arguments `min_pts` and `min_cluster_size` expose
  HDBSCAN tuning parameters that were previously hard-coded. Non-T-cell rows
  now receive `NA` (not `0`) in the `tcell_cluster_hdbscan` column.

* `scan_clustering()` — the return value is now a named list of result objects
  each containing `Lest`, `envelope`, `window_center`, and `n_cells`. New
  `nsim` and `min_cells` arguments added. Tissue-extent guard now fires before
  `seq()` is called, preventing a wrong-sign error on small datasets.

* All functions — phenotype matching now accepts both singular (`"B cell"`,
  `"T cell"`) and plural (`"B cells"`, `"T cells"`) labels via internal
  helpers `.is_bcell()` and `.is_tcell()`.

* All functions — the global `ldata` fallback now emits a `warning()` asking
  users to pass `ldata` explicitly. This behaviour is deprecated and will be
  removed in v0.3.0.

## Bug fixes

* `plot_TLS()` — replaced the Greek mu character (`\u03bc`) in axis labels
  with plain ASCII `"um"` to prevent encoding failures on systems with limited
  locale support.

* `plot_TLS()` — now returns the ggplot object truly invisibly (no implicit
  `print()` call), preventing spurious rendering during `R CMD check`.

## Package infrastructure

* Minimum R version raised to `>= 4.0.0`.
* `spatstat` imports updated to `spatstat.geom` and `spatstat.explore`.
* `VignetteBuilder: knitr` added to DESCRIPTION (was missing, causing a
  vignette WARNING in `R CMD check`).
* `LICENSE.md` and `cran-comments.md` added to `.Rbuildignore` (were
  generating a NOTE about non-standard top-level files).
* Comprehensive `testthat` unit tests covering all exported functions.
* Vignette `"tlsR-workflow"` demonstrating the full analysis pipeline.

---

# tlsR 0.1.2

* Initial CRAN release.
