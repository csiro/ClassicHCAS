test_that("reference_use returns all stages in sample-row order", {
  target <- matrix(
    c(
      0, 0, 0, 0,
      0, 0, 10, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  samples <- matrix(
    c(
      0, 0, 1, 0,
      0, 0, 2, 0,
      0, 0, 9, 0
    ),
    ncol = 4,
    byrow = TRUE
  )

  output <- reference_use(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 1,
    bin_width = 1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  expect_identical(names(output), c("id", "predicted", "density", "condition"))
  expect_equal(nrow(output), nrow(samples))
  expect_equal(output$id, 1:3)
  expect_equal(output$predicted, c(1, 2, 1))
  expect_equal(sum(output$density), 2)
  expect_equal(sum(output$condition), 2)
})

test_that("reference_use records reference-density selection", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 1,
      0, 0, 2, 0,
      0, 0, 9, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0, nrow = 20, ncol = 20)
  ref_density[2, 2] <- 0.2
  ref_density[3, 1] <- 0.9

  output <- reference_use(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 2,
    k2 = 1,
    bin_width = 1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  expect_equal(output$predicted, c(1, 1, 0))
  expect_equal(output$density, c(0, 1, 0))
  expect_equal(output$condition, c(0, 1, 0))
})

test_that("condition use optionally applies weighted-maximum attribution", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 0,
      0, 0, 2, 0
    ),
    ncol = 4,
    byrow = TRUE
  )

  output <- reference_use(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    confidence = 0.5,
    lambda = 2,
    exclude_slef = FALSE,
    num_threads = 1,
    weighted_max = TRUE,
    boost = NULL
  )

  expect_equal(output$predicted, c(1, 1))
  expect_equal(output$density, c(1, 1))
  weights <- exp(-(c(1, 2) / 2)^2)
  expected <- c(0.5, 0) + 0.5 * weights / sum(weights)
  expect_equal(output$condition, expected, tolerance = 1e-12)
  expect_equal(sum(output$condition), 1)
})

test_that("condition use attributes boosted kernel weights", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 0,
      0, 0, 2, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0, nrow = 20, ncol = 20)
  ref_density[2, 1] <- 0.8
  ref_density[3, 1] <- 1.0

  output <- reference_use(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    confidence = 0.9,
    boost = 4,
    lambda = 1,
    exclude_slef = FALSE,
    num_threads = 1,
    weighted_max = TRUE
  )

  weights <- exp(-(c(1, 2)^2))
  expected <- c(weights[1], 4 * weights[2])
  expected <- expected / sum(expected)

  expect_equal(output$condition, expected, tolerance = 1e-12)
  expect_equal(sum(output$condition), 1)
})

test_that("reference_use raster processing matches matrix processing", {
  skip_if_not_installed("terra")

  raster <- terra::rast(
    nrows = 1,
    ncols = 2,
    xmin = 0,
    xmax = 2,
    ymin = 0,
    ymax = 1,
    nlyrs = 2
  )
  terra::values(raster) <- cbind(c(0, 10), c(0, 0))
  target <- cbind(
    terra::xyFromCell(raster, 1:2),
    terra::values(raster, mat = TRUE)
  )
  samples <- matrix(
    c(
      0.5, 0.5, 1, 0,
      0.5, 0.5, 2, 0,
      0.5, 0.5, 9, 0
    ),
    ncol = 4,
    byrow = TRUE
  )
  args <- list(
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 1,
    bin_width = 1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  matrix_output <- do.call(reference_use, c(list(data = target), args))
  raster_output <- do.call(reference_use, c(list(data = raster), args))

  expect_equal(raster_output, matrix_output)
})

test_that("reference_use rejects temporal sample lists", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- list(`2000` = target, `2001` = target)

  expect_error(
    reference_use(
      data = target,
      samples = samples,
      ref_density = matrix(1, nrow = 10, ncol = 10),
      bin_width = 0.1,
      interpolate = FALSE
    ),
    "Temporal sample lists are not supported"
  )
})

test_that("reference_use aggregation is thread-consistent", {
  target <- cbind(
    x = rep(0, 20),
    y = rep(0, 20),
    predicted = seq(0, 1, length.out = 20),
    observed = seq(1, 0, length.out = 20)
  )
  samples <- cbind(
    x = rep(0, 8),
    y = rep(0, 8),
    predicted = seq(0, 1, length.out = 8),
    observed = seq(1, 0, length.out = 8)
  )
  args <- list(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 5,
    k2 = 3,
    bin_width = 0.1,
    interpolate = FALSE,
    confidence = 0.4,
    lambda = 2,
    exclude_slef = FALSE
  )

  single_thread <- do.call(reference_use, c(args, list(num_threads = 1)))
  two_threads <- do.call(reference_use, c(args, list(num_threads = 2)))

  expect_equal(two_threads$id, single_thread$id)
  expect_equal(two_threads$predicted, single_thread$predicted)
  expect_equal(two_threads$density, single_thread$density)
  expect_equal(
    two_threads$condition,
    single_thread$condition,
    tolerance = 1e-12
  )
})
