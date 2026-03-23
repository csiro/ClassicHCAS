test_that("balanced tiling exact mode can split zero-weight matrices", {
    skip_if_not_installed("terra")

    mat <- matrix(0, nrow = 4, ncol = 4)
    ext <- terra::ext(0, 4, 0, 4)

    tiles <- tiling(
        data = mat,
        n_tiles = 4,
        balanced = TRUE,
        exact = TRUE,
        spatial = TRUE,
        extent = ext
    )

    expect_s4_class(tiles, "SpatVector")
    expect_equal(nrow(tiles), 4)
})


test_that("single-row and single-column balanced tiling avoids degenerate splits", {
    skip_if_not_installed("terra")

    row_mat <- matrix(1:5, nrow = 1)
    row_ext <- terra::ext(0, 5, 0, 1)
    row_tiles <- tiling(
        data = row_mat,
        n_tiles = 5,
        balanced = TRUE,
        method = "row",
        exact = TRUE,
        spatial = TRUE,
        extent = row_ext
    )
    expect_equal(nrow(row_tiles), 5)

    col_mat <- matrix(1:5, ncol = 1)
    col_ext <- terra::ext(0, 1, 0, 5)
    col_tiles <- tiling(
        data = col_mat,
        n_tiles = 5,
        balanced = TRUE,
        method = "col",
        exact = TRUE,
        spatial = TRUE,
        extent = col_ext
    )
    expect_equal(nrow(col_tiles), 5)
})


test_that("balanced tiling covers all cells exactly once in lookup raster", {
    skip_if_not_installed("terra")

    set.seed(1)
    mat <- matrix(runif(80), nrow = 8, ncol = 10)
    ext <- terra::ext(0, 10, 0, 8)

    tiles <- tiling(
        data = mat,
        n_tiles = 7,
        balanced = TRUE,
        method = "best",
        exact = TRUE,
        spatial = TRUE,
        extent = ext
    )

    rast <- terra::rasterize(tiles, terra::rast(nrows = 8, ncols = 10, extent = ext), field = "lyr.1")
    vals <- terra::values(rast, mat = FALSE)
    expect_false(any(is.na(vals)))
    expect_equal(length(unique(vals)), 7)
})
