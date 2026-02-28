#' Number of samples within a radius
#'
#' This function calculates the number of samples (x, y coordinates) within a specified radius
#' for each pixel in a raster map.
#'
#' This function uses an integer-based distance checks for fast radius searches on either
#' geographic or projected coordinates. In geographic mode, coordinates are
#' stored in micro-degrees (degree * 1000_000) and the distance is approximated by:
#'
#'     distance² ≈ (dlat)² + (dlon × cos(lat₁))²
#'
#' where cos(lat₁) is derived from the query latitude. This avoids floating-
#' point overhead and provides substantial performance gains but introduces
#' distortion at larger distances. For applications requiring higher accuracy,
#' especially beyond regional scales (more than several 100s of kilometers in \code{radius_km}),
#' use a projected coordinate system so distances in meters can be evaluated directly.
#'
#' Ensure that \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system
#' to take advantage of parallel processing and accelerate computations. While most systems
#' include OpenMP by default, you may need to load the appropriate module if you're using an HPC
#' system.
#'
#' \strong{Note for macOS users:} Install OpenMP via Homebrew with \code{brew install libomp}
#' before installing this package.
#'
#' @param x A SpatRaster representing the study area over which sample density will be calculated.
#' @param samples_xy A matrix or data.frame containing x and y coordinates (longitude and latitude)
#' of the reference points used for density calculation.
#' @param radius_km Numeric. Specifies the search radius (buffer) in kilometers.
#' @param num_threads Integer. Specifies the number of CPU threads to be used for processing. A value
#' below 1 indicates that all available threads will be utilized. Refer to the details section for
#' more information.
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.
#'
#' @seealso \code{\link{benchmark}}
#'
#' @return A SpatRaster
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#'
#'
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
                fun = proxy_count,
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
proxy_count <- function(model, newdata, ...) {
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

