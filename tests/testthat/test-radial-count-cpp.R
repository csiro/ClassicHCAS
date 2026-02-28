test_that("radial_count returns sensible counts on a small projected raster", {
    skip_if_not_installed("terra")

    r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
    terra::values(r) <- 1

    samples_xy <- matrix(
        c(
            0.5, 0.5,
            2.5, 2.5
        ),
        ncol = 2,
        byrow = TRUE
    )

    out <- radial_count(
        x = r,
        samples_xy = samples_xy,
        radius_km = 0.0004,
        num_threads = 1L
    )

    expect_s4_class(out, "SpatRaster")

    vals <- terra::values(out, mat = FALSE)
    expect_length(vals, terra::ncell(r))
    expect_true(all(!is.na(vals)))
    expect_true(all(vals %in% c(0, 1)))
    expect_equal(sum(vals), 2)
})


test_that("radial_count is stable across thread counts", {
    skip_if_not_installed("terra")

    set.seed(42)
    r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
    terra::values(r) <- 1

    samples_xy <- cbind(
        runif(200, min = 0, max = 20),
        runif(200, min = 0, max = 20)
    )

    out_single <- radial_count(
        x = r,
        samples_xy = samples_xy,
        radius_km = 0.002,
        num_threads = 1L
    )

    out_multi <- radial_count(
        x = r,
        samples_xy = samples_xy,
        radius_km = 0.002,
        num_threads = 2L
    )

    expect_equal(
        terra::values(out_single, mat = FALSE),
        terra::values(out_multi, mat = FALSE)
    )
})
