#' Generate raster tiles using raster or matrix data
#'
#' Creates either rectangular or balanced tiles for raster or matrix data. Balanced tiles
#' attempt to distribute the data values evenly across tiles, while rectangular tiles simply
#' divide the raster into equal-sized grids. This function helps creating tiles for running
#' \code{benchmark} function over multiple tiles or systems (e.g. in a cluster). The output
#' of the \code{\link{radial_count}} function can be used to balance the run time over each tile.
#'
#' @param data A `SpatRaster` object or a numeric matrix. For balanced tiling, either a matrix
#' or raster can be used. For rectangular tiling, a raster is required.
#' @param n_tiles Integer. Number of tiles to generate.
#' @param balanced Logical. If `TRUE` (default), tiles are created to balance the sum of values within each tile.
#'   If `FALSE`, rectangular tiles of equal size are generated.
#' @param method Character. One of `"best"`, `"row"`, `"col"`, or `"both"` (default `"best"`).
#'   Specifies how the balanced tiles should be split:
#'   \describe{
#'     \item{"best"}{Automatically chooses the split direction that balances node weights and avoids overly skinny tiles.}
#'     \item{"row"}{Always splits by rows.}
#'     \item{"col"}{Always splits by columns.}
#'     \item{"both"}{Splits along both dimensions, favoring the longer dimension.}
#'   }
#' @param exact Logical. If `TRUE` (default), ensures exactly `n_tiles` are produced in balanced mode.
#' @param weighted Logical. If `TRUE` (default), tile weights are based on the scaled data values.
#'   If `FALSE`, all non-zero values are treated equally.
#' @param spatial Logical. If `TRUE`, returns a `SpatVector` polygon of tiles.
#'   If `FALSE` (default), returns a matrix or data frame with tile extents.
#' @param extent Optional. A `terra::ext` object specifying the raster extent. Required when `data` is a matrix.
#'
#' @return Either:
#' \itemize{
#'   \item A `SpatVector` of polygons representing tiles (if `spatial = TRUE`), or
#'   \item A matrix/data.frame with columns `"xmin"`, `"xmax"`, `"ymin"`, `"ymax"` for each tile (if `spatial = FALSE`).
#' }
#'
#' @details
#' - Rectangular tiling divides the raster into equal-sized tiles regardless of data values.
#' - Balanced tiling attempts to split the data such that the sum of the cell values in each tile is roughly equal.
#' - The splitting process recursively divides the matrix/raster along rows or columns, according to the `method`.
#' - NA values are treated as zero for the purpose of balanced tiling.
#'
#' @examples
#' \dontrun{
#' library(ClassicHCAS)
#' library(terra)
#'
#' r <- rast(nrows=100, ncols=100)
#' values(r) <- runif(ncell(r))
#'
#' # Balanced tiles
#' tiles_poly <- ClassicHCAS::tiling(r, n_tiles = 4, balanced = TRUE, spatial = TRUE)
#'
#' # Rectangular tiles
#' tiles_rect <- ClassicHCAS::tiling(r, n_tiles = 4, balanced = FALSE)
#'
#' # Balanced tiles from a matrix
#' mat <- as.matrix(r)
#' ext <- ext(r)
#' tiles_from_matrix <- ClassicHCAS::tiling(mat, n_tiles = 4, balanced = TRUE, extent = ext)
#' }
#'
#' @seealso
#' \code{\link[terra]{rast}}, \code{\link[terra]{as.polygons}}, \code{\link[terra]{ext}}
#'
#' @export
tiling <- function(
        data,            # matrix (only for balanced) or raster for creating the tiles
        n_tiles,         # number of tiles
        balanced = TRUE, # should it be balanced tiles, otherwise rectangular
        method = c("best", "row", "col", "both"), # balanced spliting method
        exact = TRUE,    # the exact number of tiles for balanced tiles?
        weighted = TRUE, # weights must be non-zero or NA; the NA have no weights for tiling
        spatial = FALSE, # return a polygon, otherwise a csv
        extent = NULL    # enforce raster extent to make sure it won't miss any pixels
) {
    # define the choice of splitting
    method <- match.arg(method[1], choices = c("best", "row", "col", "both"))
    n_tiles <- as.integer(n_tiles[1])
    if (is.na(n_tiles) || n_tiles < 1) stop("'n_tiles' must be a positive integer.")

    equi_tiles <- function(r, n, sp = FALSE) {
        nc <- if (n %% 2 == 1) 1 else floor(sqrt(n / 2))
        nr <- floor(n / nc)

        w <- terra::rast(terra::ext(r), nrows = nr, ncols = nc)

        if (sp) {
            return(terra::as.polygons(w))
        } else {
            return(terra::getTileExtents(r, w))
        }
    }

    if (!balanced) {
        if (any(methods::is(data, "SpatRaster"))) {
            return(
                equi_tiles(
                    r = data,
                    n = n_tiles,
                    sp = spatial
                )
            )
        } else {
            stop("Equal-sized tiles only works with rasters.")
        }
    }

    # get data and make sure there's no NA
    if (any(methods::is(data, "SpatRaster"))) {
        x <- terra::as.matrix(data, wide = TRUE)
    } else if (any(is(data, "matrix"))) {
        x <- data
    } else {
        stop("'data' must be a matrix or SpatRaster!")
    }

    # make sure every pixel is counted
    x <- ifelse(is.na(x), 0, scales::rescale(x, to = c(1, 100)))

    # make all cells weighted equally
    if (!weighted) {
        x <- ifelse(x > 0, 1, 0)
    }

    result <- tiling_cpp(
        x = x,
        n_tiles = n_tiles,
        method = method,
        exact = exact
    )

    # get the ext and make it a polygon
    if (any(methods::is(data, "SpatRaster"))) {
        if (is.null(extent)) {
            extent <- terra::ext(data)
        }
    } else {
        if (is.null(extent)) {
            stop("For matrix 'data' you need to provide the 'extent'.")
        }
    }

    outpoly <- terra::as.polygons(
        terra::rast(result, extent = extent)
    )

    if (spatial) {
        if (any(methods::is(data, "SpatRaster"))) {
            terra::crs(outpoly) <- terra::crs(data)
        }

        return(
            outpoly
        )
    } else {
        return(
            .get_tile_extent(outpoly)
        )
    }
}

# get extent of each tile
.get_tile_extent <- function(x) {
    nr <- nrow(x)
    out <- matrix(0, nrow = nr, ncol = 4)
    colnames(out) <- c("xmin", "xmax", "ymin", "ymax")
    for(i in 1:nr) {
        out[i, ] <- terra::ext(x[i, ])[1:4]
    }

    return(out)
}
