test_that("summarize_TLS returns a data.frame with expected columns", {
  data(toy_ldata, package = "tlsR")
  ldata  <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  sumtbl <- summarize_TLS(ldata)
  expect_s3_class(sumtbl, "data.frame")
  expect_true(all(c("sample", "n_TLS", "total_cells",
                    "TLS_cells", "TLS_fraction", "mean_TLS_size") %in%
                    names(sumtbl)))
})

test_that("summarize_TLS has one row per sample", {
  data(toy_ldata, package = "tlsR")
  ldata  <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  sumtbl <- summarize_TLS(ldata)
  expect_equal(nrow(sumtbl), length(ldata))
})

test_that("summarize_TLS total_cells matches nrow of sample", {
  data(toy_ldata, package = "tlsR")
  ldata  <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  sumtbl <- summarize_TLS(ldata)
  expect_equal(sumtbl$total_cells[sumtbl$sample == "ToySample"],
               nrow(ldata[["ToySample"]]))
})

test_that("summarize_TLS TLS_fraction is in [0, 1]", {
  data(toy_ldata, package = "tlsR")
  ldata  <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  sumtbl <- summarize_TLS(ldata)
  expect_true(all(sumtbl$TLS_fraction >= 0 & sumtbl$TLS_fraction <= 1,
                  na.rm = TRUE))
})

test_that("summarize_TLS warns when tls_id_knn column is absent", {
  data(toy_ldata, package = "tlsR")
  expect_warning(summarize_TLS(toy_ldata), "Run detect_TLS")
})

test_that("summarize_TLS n_TIC is populated after detect_tic", {
  data(toy_ldata, package = "tlsR")
  ldata  <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  ldata  <- detect_tic("ToySample", ldata = ldata)
  sumtbl <- summarize_TLS(ldata)
  expect_false(is.na(sumtbl$n_TIC[sumtbl$sample == "ToySample"]))
})
