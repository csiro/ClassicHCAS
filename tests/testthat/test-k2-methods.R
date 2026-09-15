.k2_method_fixture <- function() {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 1.2,
      0, 0, 2, 0.1,
      0, 0, 3, 2.9
    ),
    ncol = 4,
    byrow = TRUE
  )
  density <- matrix(0, nrow = 10, ncol = 10)
  density[2, 2] <- 0.9
  density[3, 1] <- 0.2
  density[4, 3] <- 0.1

  list(
    target = target,
    samples = samples,
    density = density,
    args = list(
      data = target,
      samples = samples,
      ref_density = density,
      radius_km = 1000,
      k1 = 3,
      k2 = 1,
      bin_width = 1,
      interpolate = FALSE,
      exclude_slef = FALSE,
      num_threads = 1
    )
  )
}


test_that("probability remains the default k2 method", {
  public_functions <- list(
    benchmark,
    reference_use,
    variable_importance,
    hcas_inspection
  )
  cpp_functions <- list(
    ClassicHCAS:::bench_cpp,
    ClassicHCAS:::reference_use_cpp,
    ClassicHCAS:::variable_importance_cpp
  )

  for (fun in c(public_functions, cpp_functions)) {
    expect_identical(formals(fun)$k2_method, "probability")
  }

  fixture <- .k2_method_fixture()
  legacy <- do.call(benchmark, fixture$args)
  explicit <- do.call(
    benchmark,
    c(fixture$args, list(k2_method = "probability"))
  )
  expect_identical(explicit, legacy)

  expect_error(
    benchmark(NULL, NULL, NULL, k2_method = "invalid"),
    "'k2_method' must be 'probability', 'residual', or 'observed'.",
    fixed = TRUE
  )
})


test_that("residual and observed methods select their requested k2 candidates", {
  fixture <- .k2_method_fixture()

  selected <- lapply(
    c("probability", "residual", "observed"),
    function(method) {
      do.call(
        reference_use,
        c(fixture$args, list(k2_method = method))
      )$density
    }
  )
  names(selected) <- c("probability", "residual", "observed")

  expect_equal(selected$probability, c(1, 0, 0))
  expect_equal(selected$residual, c(0, 0, 1))
  expect_equal(selected$observed, c(0, 1, 0))

  conditions <- vapply(
    names(selected),
    function(method) {
      unname(do.call(
        benchmark,
        c(fixture$args, list(k2_method = method))
      )[1, 1])
    },
    numeric(1)
  )
  expect_equal(conditions, c(probability = 0.9, residual = 0.1, observed = 0.2))

  signals <- vapply(
    names(selected),
    function(method) {
      unname(do.call(
        variable_importance,
        c(
          fixture$args,
          list(k2_method = method, output = "signal", boost = NULL)
        )
      )[1, 1])
    },
    numeric(1)
  )
  expect_equal(
    signals,
    c(probability = 1.2, residual = 2.9, observed = 0.1),
    tolerance = 1e-6
  )
})


test_that("k2 methods leave predicted-distance kernel weights unchanged", {
  fixture <- .k2_method_fixture()
  args <- fixture$args
  args$k2 <- 3
  args$confidence <- 0
  args$lambda <- 2
  args["boost"] <- list(NULL)

  use <- lapply(
    c("probability", "residual", "observed"),
    function(method) {
      do.call(reference_use, c(args, list(k2_method = method)))
    }
  )

  expected <- exp(-(c(1, 2, 3) / 2)^2)
  expected <- expected / sum(expected)
  for (result in use) {
    expect_equal(result$condition, expected, tolerance = 1e-12)
  }
})


test_that("new k2 methods break metric ties by sample row", {
  target <- matrix(c(0, 0, 0, 0), nrow = 1)
  samples <- matrix(
    c(
      0, 0, 1, 1.5,
      0, 0, 2, 2.5
    ),
    ncol = 4,
    byrow = TRUE
  )

  output <- reference_use(
    data = target,
    samples = samples,
    ref_density = matrix(1, nrow = 10, ncol = 10),
    radius_km = 1000,
    k1 = 2,
    k2 = 1,
    k2_method = "residual",
    bin_width = 1,
    interpolate = FALSE,
    exclude_slef = FALSE,
    num_threads = 1
  )

  expect_equal(output$density, c(1, 0))
})
