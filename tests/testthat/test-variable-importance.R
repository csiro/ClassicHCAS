# Independent R reference implementation of the per-variable metric, evaluated
# over a known selected reference set. Used to check the C++ engine.
.gaussian <- function(d, lambda) {
  exp(-(pmax(d, 0) / lambda)^2)
}

# t_obs: length-nvar observed target vector
# r_obs: n_ref x nvar matrix of observed reference values
# pd:    length-n_ref predicted distances of those references to the target
.expected_importance <- function(t_obs, r_obs, pd, lambda = 2, eps = 1e-6,
                                 boost = NULL) {
  w <- .gaussian(pd, lambda)
  if (!is.null(boost) && !is.na(boost)) {
    w[1] <- w[1] * boost
  }
  W <- sum(w)
  vapply(
    seq_len(ncol(r_obs)),
    function(v) {
      r <- r_obs[, v]
      signal <- sum(w * abs(t_obs[v] - r)) / W
      noise <- sum(outer(w, w) * abs(outer(r, r, "-"))) / W^2
      signal / (noise + eps)
    },
    numeric(1)
  )
}


test_that("single variable, equal Gaussian weights matches signal/noise", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 2,
      0, 0, 1, 4
    ),
    ncol = 4,
    byrow = TRUE
  )

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    lambda = 2,
    epsilon = 0,
    exclude_slef = FALSE,
    num_threads = 1,
    boost = NULL
  )

  # signal = (|0-2|+|0-4|)/2 = 3 ; noise = |2-4|/2 = 1 ; importance = 3
  expect_equal(dim(out), c(1L, 1L))
  expect_equal(unname(out[1, 1]), 3, tolerance = 1e-5)
})


test_that("unequal predicted distances apply Gaussian weighting", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 2,
      0, 0, 2, 4
    ),
    ncol = 4,
    byrow = TRUE
  )

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    lambda = 2,
    epsilon = 0,
    exclude_slef = FALSE,
    num_threads = 1,
    boost = NULL
  )

  expected <- .expected_importance(
    t_obs = 0,
    r_obs = matrix(c(2, 4), ncol = 1),
    pd = c(1, 2),
    lambda = 2,
    eps = 0
  )
  expect_equal(unname(out[1, 1]), expected, tolerance = 1e-5)
})


test_that("default boost weights the maximum-probability retained reference", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 2,
      0, 0, 2, 4
    ),
    ncol = 4,
    byrow = TRUE
  )
  ref_density <- matrix(0, nrow = 20, ncol = 20)
  ref_density[2, 3] <- 0.8
  ref_density[3, 5] <- 1.0

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    lambda = 1,
    epsilon = 0,
    exclude_slef = FALSE,
    num_threads = 1
  )
  unboosted <- variable_importance(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    lambda = 1,
    epsilon = 0,
    exclude_slef = FALSE,
    num_threads = 1,
    boost = NULL
  )

  expected <- .expected_importance(
    t_obs = 0,
    r_obs = matrix(c(4, 2), ncol = 1),
    pd = c(2, 1),
    lambda = 1,
    eps = 0,
    boost = 10
  )

  expect_equal(unname(out[1, 1]), expected, tolerance = 1e-5)
  expect_false(isTRUE(all.equal(out, unboosted, tolerance = 1e-5)))
})


# weighted absolute departure of target from references (the signal numerator)
.expected_signal <- function(t_obs, r_obs, pd, lambda = 2, boost = NULL) {
  w <- .gaussian(pd, lambda)
  if (!is.null(boost) && !is.na(boost)) {
    w[1] <- w[1] * boost
  }
  W <- sum(w)
  vapply(
    seq_len(ncol(r_obs)),
    function(v) sum(w * abs(t_obs[v] - r_obs[, v])) / W,
    numeric(1)
  )
}


test_that("output='signal' returns raw departure contributions and 'share' partitions them", {
  target <- matrix(c(0, 0, 0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 0.5, 0.0, 0.2, 1.0,
      0, 0, 0.0, 0.5, 0.4, 2.0,
      0, 0, 0.5, 0.5, 0.6, 5.0
    ),
    ncol = 6,
    byrow = TRUE
  )
  args <- list(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 40, ncol = 40),
    radius_km = 1000,
    k1 = 3,
    k2 = 3,
    bin_width = 0.5,
    interpolate = FALSE,
    lambda = 2,
    exclude_slef = FALSE,
    num_threads = 1,
    boost = NULL
  )

  pd <- rowSums(abs(samples[, 3:4, drop = FALSE]))
  obs <- samples[, 5:6, drop = FALSE]
  expected_signal <- .expected_signal(c(0, 0), obs, pd, 2)

  sig <- do.call(variable_importance, c(args, list(output = "signal")))
  expect_equal(as.numeric(sig[1, ]), expected_signal, tolerance = 1e-4)

  shr <- do.call(variable_importance, c(args, list(output = "share")))
  expect_equal(sum(shr[1, ]), 1, tolerance = 1e-5)
  expect_equal(
    as.numeric(shr[1, ]),
    expected_signal / sum(expected_signal),
    tolerance = 1e-4
  )
})


test_that("output='share' raster layers sum to one per cell", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 5, ncols = 5, nlyrs = 4, xmin = 0, xmax = 5,
                   ymin = 0, ymax = 5)
  set.seed(42)
  terra::values(r) <- matrix(runif(terra::ncell(r) * 4), ncol = 4)
  s <- cbind(x = runif(30, 0, 5), y = runif(30, 0, 5),
             p1 = runif(30), p2 = runif(30),
             o1 = runif(30), o2 = runif(30))

  shr <- variable_importance(r, s, matrix(1, 20, 20), radius_km = 2000,
                             k1 = 10, k2 = 5, bin_width = 0.1,
                             interpolate = FALSE, output = "share",
                             num_threads = 1)
  cell_sums <- terra::values(terra::app(shr, "sum"))
  cell_sums <- cell_sums[is.finite(cell_sums)]
  expect_true(all(abs(cell_sums - 1) < 1e-5))
})


test_that("multivariable importance matches the reference implementation", {
  # columns: x, y, pred1, pred2, obs1, obs2
  target <- matrix(c(0, 0, 0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 0.5, 0.0, 0.2, 1.0,
      0, 0, 0.0, 0.5, 0.4, 2.0,
      0, 0, 0.5, 0.5, 0.6, 5.0
    ),
    ncol = 6,
    byrow = TRUE
  )

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 40, ncol = 40),
    radius_km = 1000,
    k1 = 3,
    k2 = 3,
    bin_width = 0.5,
    interpolate = FALSE,
    lambda = 2,
    epsilon = 1e-6,
    exclude_slef = FALSE,
    num_threads = 1,
    boost = NULL
  )

  pred <- samples[, 3:4, drop = FALSE]
  obs <- samples[, 5:6, drop = FALSE]
  pd <- rowSums(abs(pred)) # target predicted is (0, 0)
  expected <- .expected_importance(
    t_obs = c(0, 0),
    r_obs = obs,
    pd = pd,
    lambda = 2,
    eps = 1e-6
  )

  expect_equal(dim(out), c(1L, 2L))
  expect_equal(as.numeric(out[1, ]), expected, tolerance = 1e-4)
})


test_that("reference-density selection drives which references contribute", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 2,
      0, 0, 1, 4,
      0, 0, 1, 6
    ),
    ncol = 4,
    byrow = TRUE
  )
  # all predicted distances are 1 (row index 2); steer selection by obs bin
  ref_density <- matrix(0, nrow = 20, ncol = 20)
  ref_density[2, 3] <- 0.9 # obs_dist 2 -> col 3
  ref_density[2, 5] <- 0.5 # obs_dist 4 -> col 5
  ref_density[2, 7] <- 0.1 # obs_dist 6 -> col 7

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 3,
    k2 = 2, # keep the two highest-probability references (obs 2 and 4)
    bin_width = 1,
    interpolate = FALSE,
    lambda = 2,
    epsilon = 0,
    exclude_slef = FALSE,
    num_threads = 1,
    boost = NULL
  )

  # kept references observed at 2 and 4 (equal weights) -> importance 3,
  # not 2.25 that including the third (obs 6) would give
  expect_equal(unname(out[1, 1]), 3, tolerance = 1e-5)
})


test_that("missing observations and empty reference sets return NaN", {
  target <- matrix(
    c(
      0, 0, 0, NA, # missing observed value
      0, 0, 0, 0   # no reference within radius
    ),
    ncol = 4,
    byrow = TRUE
  )
  samples <- matrix(c(1000, 1000, 1, 2), nrow = 1)

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1, # samples are far away -> no candidates
    k1 = 1,
    k2 = 1,
    bin_width = 1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  expect_true(is.nan(out[1, 1]))
  expect_true(is.nan(out[2, 1]))
})


test_that("raster processing matches matrix processing", {
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
  terra::values(raster) <- cbind(c(0, 0), c(0, 0))
  target <- cbind(
    terra::xyFromCell(raster, 1:2),
    terra::values(raster, mat = TRUE)
  )
  samples <- matrix(
    c(
      0.5, 0.5, 1, 2,
      0.5, 0.5, 2, 4
    ),
    ncol = 4,
    byrow = TRUE
  )
  args <- list(
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  matrix_output <- do.call(variable_importance, c(list(data = target), args))
  raster_output <- do.call(variable_importance, c(list(data = raster), args))

  expect_equal(
    unname(terra::values(raster_output)),
    unname(matrix_output),
    tolerance = 1e-5
  )
})


test_that("column names propagate from observed variables", {
  target <- matrix(
    c(0, 0, 0, 0, 0, 0),
    nrow = 1,
    dimnames = list(NULL, c("x", "y", "p_ndvi", "p_swir", "ndvi", "swir"))
  )
  samples <- matrix(
    c(
      0, 0, 0.5, 0.0, 0.2, 1.0,
      0, 0, 0.0, 0.5, 0.4, 2.0
    ),
    ncol = 6,
    byrow = TRUE,
    dimnames = list(NULL, c("x", "y", "p_ndvi", "p_swir", "ndvi", "swir"))
  )

  out <- variable_importance(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 2,
    k2 = 2,
    bin_width = 0.5,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  expect_equal(colnames(out), c("ndvi", "swir"))
})


test_that("invalid arguments are rejected", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(c(0, 0, 1, 2), nrow = 1)
  ref <- matrix(1, nrow = 20, ncol = 20)

  expect_error(
    variable_importance(target, samples, ref, k1 = 1, k2 = 2,
                        bin_width = 1, interpolate = FALSE),
    "'k2' must be less than or equal to 'k1'."
  )
  expect_error(
    variable_importance(target, samples, ref, lambda = 0,
                        bin_width = 1, interpolate = FALSE),
    "'lambda' must be one finite number greater than zero."
  )
  expect_error(
    variable_importance(target, samples, ref, epsilon = -1,
                        bin_width = 1, interpolate = FALSE),
    "'epsilon' must be one finite, non-negative number."
  )
  expect_error(
    variable_importance(target, samples = list(`2000` = target), ref,
                        bin_width = 1, interpolate = FALSE),
    "Temporal sample lists are not supported"
  )
  # benchmark()-only args must fail loudly, not be forwarded to the C++ engine
  expect_error(
    variable_importance(target, samples, ref, bin_width = 1,
                        interpolate = FALSE, confidence = 0.5),
    "does not accept benchmark\\(\\) argument"
  )
})


test_that("results are independent of thread count", {
  target <- cbind(
    x = rep(0, 25),
    y = rep(0, 25),
    p1 = seq(0, 1, length.out = 25),
    p2 = seq(1, 0, length.out = 25),
    o1 = seq(0.1, 0.9, length.out = 25),
    o2 = seq(0.9, 0.1, length.out = 25)
  )
  samples <- cbind(
    x = rep(0, 10),
    y = rep(0, 10),
    p1 = seq(0, 1, length.out = 10),
    p2 = seq(1, 0, length.out = 10),
    o1 = runif(10),
    o2 = runif(10)
  )
  args <- list(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 20, ncol = 20),
    radius_km = 1000,
    k1 = 6,
    k2 = 4,
    bin_width = 0.1,
    interpolate = FALSE,
    lambda = 2,
    exclude_slef = FALSE
  )

  one <- do.call(variable_importance, c(args, list(num_threads = 1)))
  two <- do.call(variable_importance, c(args, list(num_threads = 2)))
  expect_equal(one, two, tolerance = 1e-10)
})


test_that("aggregate_importance ranks by median and averages shares", {
  imp <- matrix(
    c(
      0.8, 0.2, 0.0,
      0.6, 0.3, 0.1,
      0.7, 0.2, 0.1
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  agg <- aggregate_importance(imp)

  expect_equal(agg$variable, c("a", "b", "c"))
  expect_equal(agg$rank, 1:3)
  expect_equal(agg$median, c(0.7, 0.2, 0.1))
  # mean of per-row shares
  shares <- imp / rowSums(imp)
  expect_equal(agg$mean_share, unname(colMeans(shares))[order(-c(0.7, 0.2, 0.1))],
               tolerance = 1e-12)
  expect_equal(unique(agg$n), 3L)
})


test_that("aggregate_importance ignores NaN cells", {
  imp <- matrix(
    c(
      0.8, 0.2,
      NaN, NaN,
      0.6, 0.4
    ),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )

  agg <- aggregate_importance(imp)
  expect_equal(agg$n[1], 2L)
  expect_equal(agg$median, c(0.7, 0.3), tolerance = 1e-12)
})
