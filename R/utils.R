# ============================================================
#  Internal helpers - not exported
# ============================================================

#' Resolve ldata: accept explicit argument or fall back to global with warning
#' @noRd
.resolve_ldata <- function(ldata) {
  if (!is.null(ldata)) {
    if (!is.list(ldata))
      stop("`ldata` must be a named list of data frames.", call. = FALSE)
    return(ldata)
  }
  if (exists("ldata", envir = .GlobalEnv, inherits = FALSE)) {
    warning(
      "Using `ldata` found in the global environment. ",
      "Passing `ldata` explicitly is strongly recommended and the global ",
      "fallback will be removed in a future version.",
      call. = FALSE
    )
    return(get("ldata", envir = .GlobalEnv, inherits = FALSE))
  }
  stop(
    "No `ldata` argument supplied and no `ldata` object found in the ",
    "global environment. Please pass your data explicitly via the `ldata` ",
    "argument.",
    call. = FALSE
  )
}

#' Validate that a sample name exists in ldata and the data frame has required
#' columns.
#' @noRd
.check_sample <- function(sample, ldata,
                           required = c("x", "y", "phenotype")) {
  if (!sample %in% names(ldata))
    stop("Sample '", sample, "' not found in `ldata`.", call. = FALSE)
  df <- ldata[[sample]]
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0)
    stop(
      "Sample '", sample, "' is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  invisible(TRUE)
}

#' Match B cells - accepts both "B cell" and "B cells" (singular and plural)
#' @noRd
.is_bcell <- function(phenotype) {
  ph <- tolower(trimws(phenotype))
  ph == "b cell" | ph == "b cells"
}

#' Match T cells - accepts both "T cell" and "T cells" (singular and plural)
#' @noRd
.is_tcell <- function(phenotype) {
  ph <- tolower(trimws(phenotype))
  ph == "t cell" | ph == "t cells"
}

#' Safe centroid of a set of (x, y) coordinates
#' @noRd
.centroid <- function(x, y) c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
