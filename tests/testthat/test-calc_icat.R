test_that("calc_icat returns a single numeric", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  val <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
  expect_length(val, 1L)
  expect_true(is.numeric(val))
})

test_that("calc_icat result is always non-negative", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  n_tls <- max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE)
  if (n_tls >= 1L) {
    vals <- vapply(seq_len(n_tls),
                   function(id) calc_icat("ToySample", tlsID = id, ldata = ldata),
                   numeric(1L))
    non_na <- vals[!is.na(vals)]
    if (length(non_na) > 0L)
      expect_true(all(non_na >= 0),
                  info = "ICAT must be non-negative for every TLS")
  } else {
    skip("No TLS detected in toy data -- non-negativity test skipped")
  }
})

test_that("calc_icat returns NA when TLS has too few cells", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  # Inject a fake TLS ID with only 2 cells
  ldata[["ToySample"]]$tls_id_knn[1:2] <- 999L
  expect_message(
    val <- calc_icat("ToySample", tlsID = 999, ldata = ldata),
    "fewer than"
  )
  expect_true(is.na(val))
})

test_that("calc_icat errors on missing tls_id_knn column", {
  data(toy_ldata, package = "tlsR")
  expect_error(
    calc_icat("ToySample", tlsID = 1, ldata = toy_ldata),
    "missing required column"
  )
})

test_that("calc_icat errors on unknown sample", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  expect_error(
    calc_icat("Ghost", tlsID = 1, ldata = ldata),
    "not found"
  )
})
