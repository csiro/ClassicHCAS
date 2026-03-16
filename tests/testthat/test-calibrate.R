test_that("calibrate clamps vector outputs to the 0 to 1 range and preserves NA", {
  out <- calibrate(
    x = c(-1, 0, 0.5, 1, 2, NA_real_),
    x_values = c(0, 1),
    y_values = c(0, 1)
  )

  expect_equal(out, c(0, 0, 0.5, 1, 1, NA_real_))
})

test_that("calibrate applies column-wise to data frame inputs", {
  x <- data.frame(
    a = c(0, 0.5, 1),
    b = c(-1, 0.5, 2)
  )

  out <- calibrate(
    x = x,
    x_values = c(0, 1),
    y_values = c(0, 1)
  )

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(3, 2))
  expect_equal(out[, 1], c(0, 0.5, 1))
  expect_equal(out[, 2], c(0, 0.5, 1))
})

test_that("calibrate validates x and y calibration vectors", {
  expect_error(
    calibrate(
      x = c(0, 1),
      x_values = "bad",
      y_values = c(0, 1)
    ),
    "x_values"
  )

  expect_error(
    calibrate(
      x = c(0, 1),
      x_values = c(0, 1),
      y_values = c(0, 0.5, 1)
    ),
    "must be the same"
  )
})
