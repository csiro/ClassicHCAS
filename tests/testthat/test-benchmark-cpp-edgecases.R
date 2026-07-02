test_that("bench_cpp handles cells with no valid neighbours without crashing", {
  raster <- matrix(c(10000, 10000, 0.1, 0.2), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 0.1, 0.2,
      10, 10, 0.15, 0.25
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0.5, nrow = 10, ncol = 10)

  output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 0.001,
    k_env = 2L,
    k_rs = 2L,
    bin_width = 0.1,
    bin_num = 10L,
    offset = 0L,
    exclude_slef = TRUE,
    make_su = FALSE,
    num_threads = 1L
  )

  expect_equal(dim(output), c(1, 1))
  expect_true(is.nan(output[1, 1]))
})

test_that("bench_cpp supports different reference density bin settings across calls", {
  raster <- matrix(c(0, 0, 0.2, 0.3), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 0.2, 0.3,
      1, 1, 0.4, 0.5
    ),
    ncol = 4,
    byrow = TRUE
  )

  ref_density_1 <- matrix(0.7, nrow = 10, ncol = 10)
  ref_density_2 <- matrix(0.2, nrow = 6, ncol = 6)

  out_1 <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density_1,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 2L,
    k_rs = 1L,
    bin_width = 0.1,
    bin_num = 10L,
    offset = 0L,
    exclude_slef = FALSE,
    make_su = FALSE,
    num_threads = 1L
  )

  out_2 <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density_2,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 2L,
    k_rs = 1L,
    bin_width = 0.2,
    bin_num = 6L,
    offset = 0L,
    exclude_slef = FALSE,
    make_su = FALSE,
    num_threads = 1L
  )

  expect_equal(out_1[1, 1], 0.7)
  expect_equal(out_2[1, 1], 0.2)
})

test_that("temporal weighting selects a year and retains its weighted probability", {
  raster <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(c(0, 0, 0.25, 0.15, 0.25), nrow = 1)
  ref_density <- matrix(0, nrow = 10, ncol = 10)
  ref_density[3, 2] <- 0.4
  ref_density[3, 3] <- 0.9

  output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 1L,
    k_rs = 1L,
    bin_width = 0.1,
    bin_num = 10L,
    offset = 0L,
    confidence = 0.5,
    exclude_slef = FALSE,
    temporal_weights = c(1, exp(-4.5)),
    make_su = FALSE,
    num_threads = 1L
  )

  expected <- 0.4
  expect_equal(output[1, 1], expected)

  weighted_output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 1L,
    k_rs = 1L,
    bin_width = 0.1,
    bin_num = 10L,
    offset = 0L,
    confidence = 0.5,
    exclude_slef = FALSE,
    temporal_weights = c(exp(-0.5), exp(-2)),
    make_su = FALSE,
    num_threads = 1L
  )

  expect_equal(weighted_output[1, 1], expected * exp(-0.5))
})

test_that("bench_cpp blends the unweighted maximum probability", {
  raster <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 0,
      0, 0, 2, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0, nrow = 10, ncol = 10)
  ref_density[2, 1] <- 0.8
  ref_density[3, 1] <- 1.0

  output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 2L,
    k_rs = 2L,
    bin_width = 1,
    bin_num = 10L,
    offset = 0L,
    confidence = 0.5,
    lambda = 2,
    exclude_slef = FALSE,
    make_su = TRUE,
    num_threads = 1L,
    boost = NULL
  )

  weights <- exp(-(c(1, 2) / 2)^2)
  probabilities <- c(0.8, 1.0)
  expected_mean <- sum(probabilities * weights) / sum(weights)
  expected_condition <-
    0.5 * max(probabilities) + 0.5 * expected_mean

  expect_equal(output[1, 1], expected_condition, tolerance = 1e-12)
  expect_equal(output[1, 2], log(sum(weights)), tolerance = 1e-12)
})

test_that("bench_cpp boosts the maximum-probability site's kernel weight", {
  raster <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 0,
      0, 0, 2, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0, nrow = 10, ncol = 10)
  ref_density[2, 1] <- 0.8
  ref_density[3, 1] <- 1.0

  output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 2L,
    k_rs = 2L,
    bin_width = 1,
    bin_num = 10L,
    offset = 0L,
    confidence = 0.9,
    lambda = 2,
    exclude_slef = FALSE,
    make_su = TRUE,
    num_threads = 1L,
    boost = 4
  )

  weights <- exp(-(c(1, 2) / 2)^2)
  probabilities <- c(0.8, 1.0)
  expected <- (
    weights[1] * probabilities[1] +
      4 * weights[2] * probabilities[2]
  ) / (weights[1] + 4 * weights[2])

  expect_equal(output[1, 1], expected, tolerance = 1e-12)
  expect_equal(output[1, 2], log(sum(weights)), tolerance = 1e-12)
})

test_that("bench_cpp stabilises weights and returns support in log space", {
  raster <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(c(0, 0, 100, 0), nrow = 1)
  ref_density <- matrix(0.4, nrow = 10, ncol = 10)

  output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 1L,
    k_rs = 1L,
    bin_width = 100,
    bin_num = 10L,
    offset = 0L,
    confidence = 0.5,
    lambda = 1,
    exclude_slef = FALSE,
    make_su = TRUE,
    num_threads = 1L
  )

  expect_equal(output[1, 1], 0.4, tolerance = 1e-12)
  expect_equal(output[1, 2], -10000, tolerance = 1e-12)
})

test_that("bench_cpp supports Cauchy weighting", {
  raster <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 0,
      0, 0, 2, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0, nrow = 10, ncol = 10)
  ref_density[2, 1] <- 0.8
  ref_density[3, 1] <- 1.0

  output <- ClassicHCAS:::bench_cpp(
    raster_vals = raster,
    sample_vals = samples,
    ref_density = ref_density,
    xy_stats = c(0, 0, 1, 1),
    radius_km = 1000,
    k_env = 2L,
    k_rs = 2L,
    bin_width = 1,
    bin_num = 10L,
    offset = 0L,
    confidence = 0.5,
    lambda = 1,
    exclude_slef = FALSE,
    make_su = TRUE,
    num_threads = 1L,
    kernel = "cauchy",
    boost = NULL
  )

  distances <- c(1, 2)
  probabilities <- c(0.8, 1.0)
  weights <- 1 / (1 + distances^2 / 1^2)
  expected_mean <- sum(probabilities * weights) / sum(weights)
  expected_condition <- 0.5 * max(probabilities) + 0.5 * expected_mean

  expect_equal(output[1, 1], expected_condition, tolerance = 1e-12)
  expect_equal(output[1, 2], log(sum(weights)), tolerance = 1e-12)
})
