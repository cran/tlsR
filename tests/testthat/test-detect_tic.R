test_that("detect_tic adds tcell_cluster_hdbscan column", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  ldata <- detect_tic("ToySample", ldata = ldata)
  expect_true("tcell_cluster_hdbscan" %in% names(ldata[["ToySample"]]))
})

test_that("detect_tic cluster IDs are non-negative integers or NA", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  ldata <- detect_tic("ToySample", ldata = ldata)
  col    <- ldata[["ToySample"]]$tcell_cluster_hdbscan
  non_na <- col[!is.na(col)]
  expect_true(all(non_na >= 0))
})

test_that("detect_tic returns ldata visibly", {
  data(toy_ldata, package = "tlsR")
  ldata  <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  result <- withVisible(detect_tic("ToySample", ldata = ldata))
  expect_true(result$visible)
})

test_that("detect_tic errors without prior detect_TLS", {
  data(toy_ldata, package = "tlsR")
  expect_error(detect_tic("ToySample", ldata = toy_ldata),
               "missing required column")
})

test_that("detect_tic messages with very large min_pts", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  expect_message(
    detect_tic("ToySample", min_pts = 99999L, ldata = ldata),
    "fewer than"
  )
})
