test_that("scan_clustering returns a list", {
  data(toy_ldata, package = "tlsR")
  result <- scan_clustering(ws = 500, sample = "ToySample",
                            phenotype = "B cells", plot = FALSE,
                            nsim = 9, ldata = toy_ldata)
  expect_type(result, "list")
})

test_that("scan_clustering list elements have expected names", {
  data(toy_ldata, package = "tlsR")
  result <- scan_clustering(ws = 500, sample = "ToySample",
                            phenotype = "B cells", plot = FALSE,
                            nsim = 9, ldata = toy_ldata)
  if (length(result) > 0) {
    elem <- result[[1]]
    expect_true(all(c("Lest", "envelope", "window_center", "n_cells") %in%
                      names(elem)))
  } else {
    skip("No significant windows found - structural test skipped")
  }
})

test_that("scan_clustering messages when tissue smaller than window", {
  # Tissue span of 3 units << ws=500 triggers the guard
  tiny <- list(Tiny = data.frame(
    x = c(1, 2, 3), y = c(1, 2, 3),
    phenotype = c("B cells", "B cells", "B cells"),
    stringsAsFactors = FALSE
  ))
  expect_message(
    scan_clustering(ws = 500, sample = "Tiny", phenotype = "B cells",
                    plot = FALSE, nsim = 9, ldata = tiny),
    "smaller than window"
  )
})

test_that("scan_clustering errors on bad phenotype argument", {
  data(toy_ldata, package = "tlsR")
  expect_error(
    scan_clustering(ws = 500, sample = "ToySample",
                    phenotype = "NK cells", plot = FALSE,
                    ldata = toy_ldata)
  )
})
