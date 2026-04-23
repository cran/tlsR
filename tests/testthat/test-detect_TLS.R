test_that("detect_TLS adds required columns", {
  data(toy_ldata, package = "tlsR")
  out <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  df  <- out[["ToySample"]]
  expect_true("tls_id_knn"   %in% names(df))
  expect_true("tls_center_x" %in% names(df))
  expect_true("tls_center_y" %in% names(df))
})

test_that("detect_TLS tls_id_knn is non-negative", {
  data(toy_ldata, package = "tlsR")
  out <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  ids <- out[["ToySample"]]$tls_id_knn
  expect_true(is.integer(ids) || is.numeric(ids))
  expect_true(all(ids >= 0, na.rm = TRUE))
})

test_that("detect_TLS detects at least one TLS in toy data", {
  data(toy_ldata, package = "tlsR")
  out <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  expect_gt(max(out[["ToySample"]]$tls_id_knn, na.rm = TRUE), 0)
})

test_that("detect_TLS errors on missing sample", {
  data(toy_ldata, package = "tlsR")
  expect_error(detect_TLS("NoSuchSample", ldata = toy_ldata), "not found")
})

test_that("detect_TLS messages when k > n B cells", {
  # Build a tiny ldata with fewer B cells than k
  tiny <- list(Tiny = data.frame(
    x = c(1, 2, 3), y = c(1, 2, 3),
    phenotype = c("B cells", "B cells", "T cells"),
    stringsAsFactors = FALSE
  ))
  expect_message(detect_TLS("Tiny", k = 30, ldata = tiny), "Not enough B cells")
})

test_that("detect_TLS respects min_B_cells filter", {
  data(toy_ldata, package = "tlsR")
  out <- detect_TLS("ToySample", k = 30, min_B_cells = 99999L,
                    ldata = toy_ldata)
  expect_equal(max(out[["ToySample"]]$tls_id_knn, na.rm = TRUE), 0)
})
