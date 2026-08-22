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
      k1 = 2,
      k2 = 1,
      bin_width = NULL,
      offset = NULL,
      interpolate = FALSE,
      num_threads = 1L
    )
  )

  expect_equal(dim(out), c(2, 1))
  expect_true(all(is.finite(out[, 1])))
})

test_that("benchmark accepts named temporal sample matrices", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  sample_2000 <- matrix(c(0, 0, 0.25, 0.15), nrow = 1)
  sample_2003 <- matrix(c(0, 0, 0.25, 0.25), nrow = 1)
  samples <- list(`2003` = sample_2003, `2000` = sample_2000)

  ref_density <- matrix(0, nrow = 10, ncol = 10)
  ref_density[3, 2] <- 0.4
  ref_density[3, 3] <- 0.9

  output <- benchmark(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 1,
    k2 = 1,
    bin_width = 0.1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    assessment_year = 2000,
    temporal_sigma = 1,
    num_threads = 1
  )

  expected <- 0.4
  expect_equal(unname(output[1, 1]), expected)

  weighted_output <- benchmark(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 1,
    k2 = 1,
    bin_width = 0.1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    assessment_year = 2001,
    temporal_sigma = 1,
    num_threads = 1
  )

  expect_equal(unname(weighted_output[1, 1]), expected * exp(-0.5))
})

test_that("temporal samples require constant predicted values", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- list(
    `2000` = matrix(c(0, 0, 0.25, 0.15), nrow = 1),
    `2001` = matrix(c(0, 0, 0.30, 0.15), nrow = 1)
  )

  expect_error(
    benchmark(
      data = target,
      samples = samples,
      ref_density = matrix(1, nrow = 10, ncol = 10),
      bin_width = 0.1,
      interpolate = FALSE,
      assessment_year = 2000,
      temporal_sigma = 1
    ),
    "Predicted RS values must be constant"
  )
})

test_that("benchmark interpolates reference_density objects", {
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

  expect_no_error(
    out <- benchmark(
      data = samples[1:2, , drop = FALSE],
      samples = samples,
      ref_density = ref_norm,
      radius_km = 1000,
      k1 = 2,
      k2 = 1,
      bin_width = NULL,
      offset = NULL,
      interpolate = TRUE,
      num_threads = 1L
    )
  )

  expect_equal(dim(out), c(2, 1))
  expect_true(all(is.finite(out[, 1])))
})
test_that("public filter arguments are named k1 and k2", {
  public_functions <- list(benchmark, reference_use, variable_importance)

  for (fun in public_functions) {
    arguments <- names(formals(fun))
    expect_true(all(c("k1", "k2") %in% arguments))
    expect_false(any(c("k_pred", "k_obs") %in% arguments))
  }

  expect_error(
    benchmark(NULL, NULL, NULL, k_pred = 1),
    "were renamed to 'k1' and 'k2'",
    fixed = TRUE
  )
  expect_error(
    variable_importance(NULL, NULL, NULL, k_obs = 1),
    "were renamed to 'k1' and 'k2'",
    fixed = TRUE
  )
})

test_that("kernel bandwidth defaults to lambda one", {
  public_functions <- list(
    benchmark,
    reference_use,
    variable_importance,
    hcas_inspection
  )

  for (fun in public_functions) {
    expect_identical(formals(fun)$lambda, 1.0)
  }

  expect_identical(formals(ClassicHCAS:::bench_cpp)$lambda, 1.0)
  expect_identical(formals(ClassicHCAS:::reference_use_cpp)$lambda, 1.0)
  expect_identical(formals(ClassicHCAS:::variable_importance_cpp)$lambda, 1.0)
})

test_that("benchmark uses the unweighted maximum without a public flag", {
  expect_false("weighted_max" %in% names(formals(benchmark)))
  expect_false("weighted_max" %in% names(formals(ClassicHCAS:::bench_cpp)))
  expect_identical(formals(reference_use)$weighted_max, FALSE)
  expect_identical(formals(ClassicHCAS:::reference_use_cpp)$weighted_max, FALSE)
  expect_false("weighted_max" %in% names(formals(hcas_inspection)))

  expect_error(
    benchmark(NULL, NULL, NULL, weighted_max = NA),
    "'weighted_max' is not an argument to benchmark().",
    fixed = TRUE
  )
})

test_that("public boost defaults require a positive finite factor", {
  public_functions <- list(
    benchmark,
    reference_use,
    variable_importance,
    hcas_inspection
  )
  for (fun in public_functions) {
    expect_identical(formals(fun)$k1, 70)
    expect_identical(formals(fun)$k2, 10)
    expect_identical(formals(fun)$boost, quote(k2))
  }
  expect_identical(formals(ClassicHCAS:::bench_cpp)$k_env, 70L)
  expect_identical(formals(ClassicHCAS:::bench_cpp)$k_rs, 10L)
  expect_identical(formals(ClassicHCAS:::reference_use_cpp)$k_env, 70L)
  expect_identical(formals(ClassicHCAS:::reference_use_cpp)$k_rs, 10L)
  expect_identical(formals(ClassicHCAS:::variable_importance_cpp)$k_env, 70L)
  expect_identical(formals(ClassicHCAS:::variable_importance_cpp)$k_rs, 10L)
  expect_true("boost" %in% names(formals(ClassicHCAS:::bench_cpp)))
  expect_true("boost" %in% names(formals(ClassicHCAS:::reference_use_cpp)))
  expect_true("boost" %in% names(formals(ClassicHCAS:::variable_importance_cpp)))
  expect_null(ClassicHCAS:::.check_boost(NULL))
  expect_null(ClassicHCAS:::.check_boost(NA_real_))

  expect_error(
    benchmark(NULL, NULL, NULL, boost = 0),
    "'boost' must be NULL, NA, or one finite number greater than zero.",
    fixed = TRUE
  )
  expect_error(
    reference_use(NULL, NULL, NULL, boost = Inf),
    "'boost' must be NULL, NA, or one finite number greater than zero.",
    fixed = TRUE
  )
  expect_error(
    variable_importance(NULL, NULL, NULL, boost = 0),
    "'boost' must be NULL, NA, or one finite number greater than zero.",
    fixed = TRUE
  )
})

test_that("Gaussian kernel is the default", {
  public_functions <- list(
    benchmark,
    reference_use,
    variable_importance,
    hcas_inspection
  )
  for (fun in public_functions) {
    expect_identical(
      eval(formals(fun)$kernel),
      c("Gaussian", "Cauchy")
    )
  }

  expect_identical(formals(ClassicHCAS:::bench_cpp)$kernel, "gaussian")
  expect_identical(
    formals(ClassicHCAS:::reference_use_cpp)$kernel,
    "gaussian"
  )
  expect_identical(
    formals(ClassicHCAS:::variable_importance_cpp)$kernel,
    "gaussian"
  )
  expect_identical(ClassicHCAS:::.check_kernel("Gaussian"), "gaussian")
  expect_identical(ClassicHCAS:::.check_kernel("cauchy"), "cauchy")

  expect_error(
    benchmark(NULL, NULL, NULL, kernel = "invalid"),
    "'kernel' must be 'Gaussian'/'gaussian' or 'Cauchy'/'cauchy'.",
    fixed = TRUE
  )
})

test_that("temporal_sigma controls the experimental temporal mode", {
  arguments <- formals(benchmark)
  expect_false("temporal_correct" %in% names(arguments))
  expect_null(arguments$temporal_sigma)
  expect_false("temporal_weighted" %in% names(arguments))
  expect_false("temporal_weighted" %in% names(formals(ClassicHCAS:::bench_cpp)))

  target <- matrix(c(0, 0, 0.25, 0.15), nrow = 1)
  samples <- matrix(c(0, 0, 0.25, 0.15), nrow = 1)
  ref_density <- matrix(1, nrow = 10, ncol = 10)

  default_output <- benchmark(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 1,
    k2 = 1,
    bin_width = 0.1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )
  na_output <- benchmark(
    data = target,
    samples = samples,
    ref_density = ref_density,
    radius_km = 1000,
    k1 = 1,
    k2 = 1,
    bin_width = 0.1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    temporal_sigma = NA,
    num_threads = 1
  )

  expect_equal(na_output, default_output)
  expect_error(
    benchmark(NULL, NULL, NULL, temporal_correct = TRUE),
    "'temporal_correct' has been removed",
    fixed = TRUE
  )
  expect_error(
    benchmark(NULL, NULL, NULL, temporal_weighted = NA),
    "'temporal_weighted' has been removed",
    fixed = TRUE
  )
})
