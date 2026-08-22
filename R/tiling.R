#' Create raster processing tiles
#'
#' Creates rectangular or value-balanced tiles for splitting large raster
#' processing jobs, especially HCAS benchmarking runs, into smaller spatial
#' chunks.
#'
#' @details
#' Large HCAS benchmarking jobs can be uneven because areas with many nearby
#' reference samples require more work than sparse areas. \code{tiling()} can
#' split a raster or matrix into tiles whose total cell weights are roughly
#' balanced, making distributed or multi-node processing more even. A common
#' workflow is to run \code{\link{radial_count}} first, then use the resulting
#' count raster as the \code{data} argument for balanced tiling.
#'
#' If \code{balanced = FALSE}, the function creates simple rectangular tiles of
#' similar size and ignores cell values. Rectangular tiling requires a
#' \pkg{terra} \code{SpatRaster}. If \code{balanced = TRUE}, the function treats
#' \code{NA} as zero, rescales non-missing values to positive weights, and
#' recursively splits the raster or matrix so each tile has a similar total
#' weight. If \code{weighted = FALSE}, all non-zero cells are treated equally.
#'
#' @param data A \pkg{terra} \code{SpatRaster} or numeric matrix. Balanced tiling
#' accepts either form; rectangular tiling requires a raster.
#' @param n_tiles Integer. Number of tiles to generate.
#' @param balanced Logical. If \code{TRUE}, create tiles with approximately
#' balanced total weights. If \code{FALSE}, create rectangular tiles of similar
#' size.
#' @param method Character. Splitting strategy for balanced tiles. One of
#' \code{"best"}, \code{"row"}, \code{"col"}, or \code{"both"}:
#' \describe{
#'   \item{\code{"best"}}{Automatically chooses a split direction that balances
#'   weights while avoiding very narrow tiles.}
#'   \item{\code{"row"}}{Always splits by rows.}
#'   \item{\code{"col"}}{Always splits by columns.}
#'   \item{\code{"both"}}{Splits along both dimensions, favouring the longer
#'   dimension.}
#' }
#' @param exact Logical. If \code{TRUE}, force exactly \code{n_tiles} in balanced
#' mode.
#' @param weighted Logical. If \code{TRUE}, tile weights are based on the scaled
#' data values. If \code{FALSE}, all non-zero cells receive equal weight.
#' @param spatial Logical. If \code{TRUE}, return a \pkg{terra}
#' \code{SpatVector} polygon layer. If \code{FALSE}, return tile extents.
#' @param extent Optional \code{\link[terra]{ext}} object specifying the raster
#' extent. Required when \code{data} is a matrix.
#'
#' @return If \code{spatial = TRUE}, a \pkg{terra} \code{SpatVector} of tile
#' polygons. Otherwise, a matrix with columns \code{xmin}, \code{xmax},
#' \code{ymin}, and \code{ymax}.
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' r <- terra::rast(nrows = 20, ncols = 20)
#' terra::values(r) <- runif(terra::ncell(r))
#'
#' # Balanced tiles from raster weights.
#' balanced_tiles <- tiling(r, n_tiles = 4, balanced = TRUE)
#'
#' # Rectangular tiles.
#' rectangular_tiles <- tiling(r, n_tiles = 4, balanced = FALSE)
#'
#' # Balanced tiles from a matrix need an explicit extent.
#' mat <- terra::as.matrix(r, wide = TRUE)
#' matrix_tiles <- tiling(mat, n_tiles = 4, extent = terra::ext(r))
#' }
#'
#' @seealso
#' \code{\link{radial_count}}, \code{\link{benchmark}},
#' \code{\link[terra]{rast}}, \code{\link[terra]{as.polygons}},
#' \code{\link[terra]{ext}}
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
    } else if (is.matrix(data)) {
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
