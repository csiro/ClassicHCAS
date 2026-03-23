test_that("ref_density drop_features matches manual feature removal", {
  samples_full <- matrix(
    c(
      0, 0, 0.20, 0.30, 0.22, 0.31,
      1, 1, 0.35, 0.45, 0.37, 0.46,
      2, 2, 0.50, 0.60, 0.52, 0.62,
      3, 3, 0.65, 0.75, 0.67, 0.78
    ),
    ncol = 6,
    byrow = TRUE
  )
  samples_reduced <- samples_full[, c(1, 2, 4, 6), drop = FALSE]

  dropped <- ref_density(
    data = samples_full,
    radius_km = 1000,
    bin_width = 0.1,
    bin_num = 10,
    drop_features = 1,
    num_threads = 1L
  )
  reduced <- ref_density(
    data = samples_reduced,
    radius_km = 1000,
    bin_width = 0.1,
    bin_num = 10,
    num_threads = 1L
  )

  expect_equal(unclass(dropped), unclass(reduced))
  expect_equal(attr(dropped, "bin.width"), attr(reduced, "bin.width"))
})

test_that("benchmark drop_features matches manual feature removal", {
  samples_full <- matrix(
    c(
      0, 0, 0.20, 0.30, 0.21, 0.31,
      1, 1, 0.40, 0.55, 0.41, 0.56,
      2, 2, 0.60, 0.70, 0.61, 0.72,
      3, 3, 0.80, 0.90, 0.81, 0.91
    ),
    ncol = 6,
    byrow = TRUE
  )
  samples_reduced <- samples_full[, c(1, 2, 4, 6), drop = FALSE]
  targets_full <- samples_full[1:2, , drop = FALSE]
  targets_reduced <- samples_reduced[1:2, , drop = FALSE]

  ref_drop <- ref_density(
    data = samples_full,
    radius_km = 1000,
    bin_width = 0.1,
    bin_num = 10,
    drop_features = 1,
    num_threads = 1L
  )
  ref_reduced <- ref_density(
    data = samples_reduced,
    radius_km = 1000,
    bin_width = 0.1,
    bin_num = 10,
    num_threads = 1L
  )

  dropped <- benchmark(
    data = targets_full,
    samples = samples_full,
    ref_density = unclass(ref_drop),
    radius_km = 1000,
    k_pred = 2,
    k_obs = 1,
    bin_width = 0.1,
    offset = 0,
    interpolate = FALSE,
    drop_features = 1,
    num_threads = 1L
  )
  reduced <- benchmark(
    data = targets_reduced,
    samples = samples_reduced,
    ref_density = unclass(ref_reduced),
    radius_km = 1000,
    k_pred = 2,
    k_obs = 1,
    bin_width = 0.1,
    offset = 0,
    interpolate = FALSE,
    num_threads = 1L
  )

  expect_equal(dropped, reduced)
})

test_that("drop_features uses 1-based RS feature positions", {
  samples <- matrix(
    c(
      0, 0, 0.20, 0.30, 0.22, 0.31,
      1, 1, 0.35, 0.45, 0.37, 0.46,
      2, 2, 0.50, 0.60, 0.52, 0.62,
      3, 3, 0.65, 0.75, 0.67, 0.78
    ),
    ncol = 6,
    byrow = TRUE
  )

  expect_error(
    ref_density(
      data = samples,
      radius_km = 1000,
      bin_width = 0.1,
      bin_num = 10,
      drop_features = 3,
      num_threads = 1L
    ),
    "'drop_features' must be between 1 and 2\\."
  )
})
