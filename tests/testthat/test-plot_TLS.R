test_that("plot_TLS returns a ggplot object invisibly", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  # Use withVisible to capture return value without rendering to screen
  result <- withVisible(plot_TLS("ToySample", ldata = ldata))
  expect_s3_class(result$value, "ggplot")
  expect_false(result$visible)  # must be invisible
})

test_that("plot_TLS returns ggplot after detect_tic", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  ldata <- detect_tic("ToySample", ldata = ldata)
  result <- withVisible(plot_TLS("ToySample", ldata = ldata, show_tic = TRUE))
  expect_s3_class(result$value, "ggplot")
})

test_that("plot_TLS errors on missing sample", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  expect_error(plot_TLS("Ghost", ldata = ldata), "not found")
})
