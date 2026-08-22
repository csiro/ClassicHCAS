test_that("hcas_inspection builds a Shiny application", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyWidgets")
  skip_if_not_installed("leaflet")
  skip_if_not_installed("ggplot2")

  samples <- matrix(
    c(
      150.0, -35.0, 0.1, 0.1,
      150.1, -35.0, 0.2, 0.2,
      150.2, -35.0, 0.3, 0.3
    ),
    ncol = 4,
    byrow = TRUE
  )
  density <- matrix(1, nrow = 10, ncol = 10)
  class(density) <- c("reference_density", "matrix", "array")
  attr(density, "bin.width") <- 0.1
  attr(density, "offset") <- 0L

  app <- hcas_inspection(
    data = samples,
    samples = samples,
    ref_density = density,
    radius_km = 100,
    k1 = 2,
    k2 = 1,
    interpolate = FALSE,
    launch = FALSE
  )

  expect_s3_class(app, "shiny.appobj")
})

test_that("hcas_inspection builds with raster data and xy-only samples", {
  skip_if_not_installed("terra")
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyWidgets")
  skip_if_not_installed("leaflet")
  skip_if_not_installed("ggplot2")

  raster <- terra::rast(
    nrows = 1,
    ncols = 3,
    xmin = 150,
    xmax = 150.3,
    ymin = -35.1,
    ymax = -35.0,
    nlyrs = 2,
    crs = "EPSG:4326"
  )
  terra::values(raster) <- cbind(c(0.1, 0.2, 0.3), c(0.1, 0.2, 0.3))
  samples_xy <- terra::xyFromCell(raster, 1:3)
  density <- matrix(1, nrow = 10, ncol = 10)
  class(density) <- c("reference_density", "matrix", "array")
  attr(density, "bin.width") <- 0.1
  attr(density, "offset") <- 0L

  app <- hcas_inspection(
    data = raster,
    samples = samples_xy,
    ref_density = density,
    radius_km = 100,
    k1 = 2,
    k2 = 1,
    interpolate = FALSE,
    launch = FALSE
  )

  expect_s3_class(app, "shiny.appobj")
})

test_that("hcas_inspection uses explicit CRS over non-transformable raster CRS", {
  skip_if_not_installed("terra")
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyWidgets")
  skip_if_not_installed("leaflet")
  skip_if_not_installed("ggplot2")

  raster <- terra::rast(
    nrows = 1,
    ncols = 1,
    xmin = -1920015,
    xmax = -1919565,
    ymin = -4896135,
    ymax = -4895685,
    nlyrs = 2
  )
  terra::values(raster) <- matrix(c(0.1, 0.1), nrow = 1)
  terra::crs(raster) <- paste0(
    "ENGCRS[\"GDA94 / Australian Albers\",",
    "EDATUM[\"Unknown engineering datum\"],",
    "CS[Cartesian,2],",
    "AXIS[\"Easting (E)\",east,ORDER[1],LENGTHUNIT[\"metre\",1]],",
    "AXIS[\"Northing (N)\",north,ORDER[2],LENGTHUNIT[\"metre\",1]]]"
  )
  samples <- cbind(
    terra::xyFromCell(raster, 1),
    matrix(c(0.1, 0.1), nrow = 1)
  )
  density <- matrix(1, nrow = 10, ncol = 10)
  class(density) <- c("reference_density", "matrix", "array")
  attr(density, "bin.width") <- 0.1
  attr(density, "offset") <- 0L

  app <- hcas_inspection(
    data = raster,
    samples = samples,
    ref_density = density,
    radius_km = 100,
    k1 = 1,
    k2 = 1,
    interpolate = FALSE,
    crs = "EPSG:3577",
    launch = FALSE
  )

  expect_s3_class(app, "shiny.appobj")
})

test_that("inspection raster target snaps coordinates to cell centre", {
  skip_if_not_installed("terra")

  raster <- terra::rast(
    nrows = 2,
    ncols = 2,
    xmin = 0,
    xmax = 2,
    ymin = 0,
    ymax = 2,
    nlyrs = 2,
    crs = "EPSG:4326"
  )
  terra::values(raster) <- cbind(1:4, 5:8)

  target <- ClassicHCAS:::.inspection_raster_target(raster, 1.8, 0.2)

  expect_equal(
    unname(target[1, 1:2]),
    unname(as.numeric(terra::xyFromCell(raster, 4)))
  )
  expect_equal(unname(target[1, 3:4]), c(4, 8))
})

test_that("inspection raster default target is the extent centre", {
  skip_if_not_installed("terra")

  raster <- terra::rast(
    nrows = 2,
    ncols = 2,
    xmin = 0,
    xmax = 2,
    ymin = 0,
    ymax = 2,
    nlyrs = 2,
    crs = "EPSG:4326"
  )
  # The centre need not fall on a populated cell; every cell here is empty.
  terra::values(raster) <- NA

  xy <- ClassicHCAS:::.inspection_raster_default_xy(raster)

  expect_equal(unname(xy), c(1, 1))
})

test_that("inspection point evaluates an extracted raster target", {
  skip_if_not_installed("terra")

  raster <- terra::rast(
    nrows = 1,
    ncols = 3,
    xmin = 150,
    xmax = 150.3,
    ymin = -35.1,
    ymax = -35.0,
    nlyrs = 2,
    crs = "EPSG:4326"
  )
  terra::values(raster) <- cbind(c(0.1, 0.2, 0.3), c(0.1, 0.2, 0.3))
  samples <- cbind(
    terra::xyFromCell(raster, 1:3),
    terra::values(raster, mat = TRUE)
  )
  target <- ClassicHCAS:::.inspection_raster_target(raster, 150.21, -35.05)

  result <- ClassicHCAS:::.inspection_point(
    target = target,
    samples = samples,
    ref_density = matrix(1, nrow = 10, ncol = 10),
    xy_stats = c(0, 0, 1, 1),
    xy_penalty = 0,
    radius_km = 100,
    k1 = 2,
    k2 = 1,
    bin_width = 0.1,
    offset = 0L,
    confidence = 0.5,
    boost = 10,
    lambda = 1,
    exclude_slef = FALSE,
    drop_features = NULL,
    num_threads = 1,
    kernel = "gaussian",
    geographic = TRUE,
    crs = NULL
  )

  expect_true(is.finite(result$condition))
  expect_s3_class(result$nearby, "data.frame")
  expect_s3_class(result$selected, "data.frame")
  expect_equal(
    unname(unlist(result$target[1, c("x", "y")])),
    unname(target[1, 1:2])
  )
})

test_that("inspection raster click coordinates snap to raster cell centre", {
  skip_if_not_installed("terra")

  raster <- terra::rast(
    nrows = 2,
    ncols = 2,
    xmin = 145,
    xmax = 147,
    ymin = -40,
    ymax = -38,
    nlyrs = 2,
    crs = "EPSG:4326"
  )
  terra::values(raster) <- 1

  xy <- ClassicHCAS:::.inspection_raster_click_xy(
    raster,
    lng = 146.8,
    lat = -39.8,
    crs = terra::crs(raster)
  )

  expect_equal(unname(xy), unname(as.numeric(terra::xyFromCell(raster, 4))))
})
