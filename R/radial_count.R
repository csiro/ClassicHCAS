#' Count samples within a radius around each raster cell
#'
#' Counts how many sample points fall within a specified radius of each cell in a
#' raster. In HCAS workflows this is mainly an operational helper for mapping
#' local reference-sample support and for creating workload weights before
#' tiling large benchmarking jobs.
#'
#' @details
#' \code{radial_count()} does not calculate habitat condition. It produces a
#' sample-density raster that can help diagnose sparse reference coverage or
#' guide \code{\link{tiling}} so densely sampled areas with heavy computation
#' are balanced across tiles.
#'
#' In geographic coordinates, radius searches use a fast integer approximation.
#' Coordinates are stored in micro-degrees (\code{degree * 1000000}) and distance
#' is approximated by:
#'
#' \deqn{distance^2 \approx dlat^2 + (dlon \times \cos(lat_1))^2}
#'
#' where \eqn{\cos(lat_1)} is derived from the query latitude. This is efficient
#' for large analyses but introduces distortion over large areas. For high
#' accuracy at broad regional or continental radii, use a projected coordinate
#' reference system so distances can be evaluated in metres.
#'
#' \code{num_threads} uses OpenMP when available. On macOS, installing OpenMP
#' support with \code{brew install libomp} before installing the package may be
#' required for multi-threaded execution.
#'
#' @param x A \pkg{terra} \code{SpatRaster} whose cells define the locations
#' where sample counts are calculated. Only the first layer is used.
#' @param samples_xy A two-column matrix or data.frame containing sample
#' coordinates in the same coordinate reference system as \code{x}.
#' @param radius_km Numeric. Search radius, in kilometres.
#' @param num_threads Integer. Number of CPU threads to use. Values below 1 use
#' all available OpenMP threads.
#' @param ... Additional arguments passed to \code{\link[terra]{interpolate}},
#' such as \code{filename}, \code{overwrite}, or \code{wopt}.
#'
#' @seealso \code{\link{benchmark}}, \code{\link{tiling}}
#'
#' @return A \pkg{terra} \code{SpatRaster} containing sample counts.
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' r <- terra::rast(
#'     nrows = 10, ncols = 10,
#'     xmin = 0, xmax = 1, ymin = 0, ymax = 1,
#'     crs = "EPSG:4326"
#' )
#'
#' samples <- cbind(x = c(0.2, 0.8), y = c(0.2, 0.8))
#' counts <- radial_count(r, samples, radius_km = 50, num_threads = 1)
#' counts
#' }
radial_count <- function(
        x,
        samples_xy,
        radius_km = 200,
        num_threads = -1,
        ...) {

    # check samples
    if (.is_mat(samples_xy)) {
        samples_xy <- .check_mat(samples_xy)
    } else {
        stop("'samples_xy' must be a matrix or an object convertibe to matrix.")
    }

    # check terra is available
    .check_pkgs("terra")
    # get the raster layers
    x <- .check_rast(x)

    tryCatch(
        {
            output <- terra::interpolate(
                object = x[[1]],
                model = list(),
                fun = .proxy_count,
                xy = samples_xy[, 1:2],
                radius_km = radius_km,
                geographic = .is_lonlat(x),
                num_threads = num_threads,
                ...
            )
        },
        error = function(cond) {
            stop("Radial count calculation failed!\n", cond)
        }
    )

    return(
        output
    )
}

# wrapper function for radial_count_cpp
.proxy_count <- function(model, newdata, ...) {
    nr <- nrow(newdata)

    tryCatch(
        {
            pcount <- radial_count_cpp(
                rast = as.matrix(newdata),
                ...
            )
        },
        error = function(cond) {
            message("Error: the radial_count C++ function faild, returning -2!")
            # return error values -0.02
            return(
                rep(-2, nr)
            )
        }
    )

    return(pcount)
}

