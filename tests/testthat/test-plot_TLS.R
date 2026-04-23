test_that("plot_TLS returns a ggplot object invisibly", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  result <- withVisible(plot_TLS("ToySample", ldata = ldata))
  expect_s3_class(result$value, "ggplot")
  expect_false(result$visible)
})

test_that("plot_TLS accepts bg_alpha and tic_size_mult arguments", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  # Should not error with non-default values
  expect_silent(
    plot_TLS("ToySample", ldata = ldata,
             bg_alpha = 0.1, tic_size_mult = 2.5)
  )
})

test_that("plot_TLS returns ggplot after detect_tic", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  ldata <- detect_tic("ToySample", ldata = ldata)
  result <- withVisible(
    plot_TLS("ToySample", ldata = ldata, show_tic = TRUE)
  )
  expect_s3_class(result$value, "ggplot")
})

test_that("plot_TLS bg_alpha is lower than alpha by default", {
  # Verify the documented defaults hold
  formals_list <- formals(plot_TLS)
  expect_lt(formals_list$bg_alpha, formals_list$alpha)
})

test_that("plot_TLS errors on missing sample", {
  data(toy_ldata, package = "tlsR")
  ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
  expect_error(plot_TLS("Ghost", ldata = ldata), "not found")
})
