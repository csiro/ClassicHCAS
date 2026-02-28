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
