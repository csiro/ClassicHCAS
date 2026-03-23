test_that("benchmark uses reference density attributes when they are not supplied", {
  skip_if_not_installed("terra")

  samples <- matrix(
    c(
      0, 0, 0.2, 0.3, 0.2, 0.3,
      1, 1, 0.4, 0.5, 0.4, 0.5,
      2, 2, 0.6, 0.7, 0.6, 0.7,
      3, 3, 0.8, 0.9, 0.8, 0.9
    ),
    ncol = 6,
    byrow = TRUE
  )

  ref <- ref_density(
    data = samples,
    radius_km = 1000,
    bin_width = 0.1,
    bin_num = 10,
    num_threads = 1L
  )
  ref_norm <- normalise(
    x = ref,
    trim_size = 6,
    offset = 1
  )

  expect_no_warning(
    out <- benchmark(
      data = samples[1:2, , drop = FALSE],
      samples = samples,
      ref_density = ref_norm,
      radius_km = 1000,
      k_pred = 2,
      k_obs = 1,
      bin_width = NULL,
      offset = NULL,
      interpolate = FALSE,
      num_threads = 1L
    )
  )

  expect_equal(dim(out), c(2, 1))
  expect_true(all(is.finite(out[, 1])))
})
