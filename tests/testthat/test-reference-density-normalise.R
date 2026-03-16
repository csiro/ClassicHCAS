test_that("ref_density supports a single remote-sensing feature", {
  skip_if_not_installed("terra")

  data <- matrix(
    c(
      0, 0, 0.2, 0.2,
      1, 1, 0.4, 0.4,
      2, 2, 0.6, 0.6,
      3, 3, 0.8, 0.8
    ),
    ncol = 4,
    byrow = TRUE
  )

  out <- ref_density(
    data = data,
    radius_km = 1000,
    bin_width = 0.1,
    bin_num = 10,
    num_threads = 1L
  )

  expect_true(inherits(out, "reference_density"))
  expect_equal(dim(out), c(10, 10))
  expect_equal(attr(out, "bin.width"), 0.1)
})

test_that("normalise handles reference_density objects and normalises rows", {
  skip_if_not_installed("terra")

  x <- matrix(seq_len(36), nrow = 6)
  class(x) <- c("reference_density", "matrix", "array")
  attr(x, "bin.width") <- 0.5

  outfile <- tempfile("normalised-reference-density-")
  out <- normalise(
    x = x,
    trim_size = 4,
    offset = -3,
    filename = outfile
  )

  expect_true(inherits(out, "reference_density"))
  expect_equal(dim(out), c(4, 4))
  expect_equal(attr(out, "bin.width"), 0.5)
  expect_equal(attr(out, "offset"), 0)
  expect_equal(rowSums(out), rep(1, 4), tolerance = 1e-8)
  expect_true(file.exists(paste0(outfile, ".txt")))
})
