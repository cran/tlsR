test_that("scan_clustering returns a list", {
  data(toy_ldata, package = "tlsR")
  result <- scan_clustering(ws = 500, sample = "ToySample",
                            phenotype = "B cells", plot = FALSE,
                            ldata = toy_ldata)
  expect_type(result, "list")
})

test_that("scan_clustering B-cell results are Lest objects", {
  data(toy_ldata, package = "tlsR")
  result <- scan_clustering(ws = 200, sample = "ToySample",
                            phenotype = "B cells", plot = FALSE,
                            ldata = toy_ldata)
  if (length(result$B) > 0) {
    expect_true(inherits(result$B[[1]], "fv"))
  } else {
    skip("No B-cell windows analysed -- structural test skipped")
  }
})

test_that("scan_clustering T-cell results are Lest objects", {
  data(toy_ldata, package = "tlsR")
  result <- scan_clustering(ws = 200, sample = "ToySample",
                            phenotype = "T cells", plot = FALSE,
                            ldata = toy_ldata)
  if (length(result$T) > 0) {
    expect_true(inherits(result$T[[1]], "fv"))
  } else {
    skip("No T-cell windows analysed -- structural test skipped")
  }
})

test_that("scan_clustering Both returns B and T elements", {
  data(toy_ldata, package = "tlsR")
  result <- scan_clustering(ws = 200, sample = "ToySample",
                            phenotype = "Both", plot = FALSE,
                            ldata = toy_ldata)
  # At least one of B or T should be present in toy data
  expect_true("B" %in% names(result) || "T" %in% names(result))
})

test_that("scan_clustering: nsim argument is no longer accepted", {
  data(toy_ldata, package = "tlsR")
  expect_error(
    scan_clustering(ws = 500, sample = "ToySample",
                    phenotype = "B cells", plot = FALSE,
                    nsim = 9, ldata = toy_ldata)
  )
})

test_that("scan_clustering messages when no windows meet thresholds", {
  data(toy_ldata, package = "tlsR")
  expect_message(
    scan_clustering(ws = 500, sample = "ToySample",
                    phenotype = "B cells", plot = FALSE,
                    min_cells = 1e6L, ldata = toy_ldata),
    "no windows"
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

test_that("scan_clustering returns empty list when no phenotype cells present", {
  no_b <- list(NoBcells = data.frame(
    x = runif(100, 0, 2000), y = runif(100, 0, 2000),
    phenotype = rep("T cells", 100),
    stringsAsFactors = FALSE
  ))
  expect_message(
    result <- scan_clustering(ws = 500, sample = "NoBcells",
                              phenotype = "B cells", plot = FALSE,
                              ldata = no_b),
    "no windows"
  )
  expect_equal(length(result), 0L)
})
