#' Number of samples within a proximity (radial search)
#'
#' This function calculates the number of samples (x, y coordinates) within a specified radius
#' for each pixel in a raster map, using Euclidean distance.
#'
#' Ensure that \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system
#' to take advantage of parallel processing and accelerate computations. While most systems
#' include OpenMP by default, you may need to load the appropriate module if you're using an HPC
#' system.
#'
#' @param x A SpatRaster representing the study area over which sample density will be calculated.
#' @param samples_xy A matrix or data.frame containing x and y coordinates (longitude and latitude)
#' of the reference points used for density calculation.
#' @param radius_km Numeric. Specifies the search radius (buffer) in kilometers.
#' @param num_threads Integer. Specifies the number of CPU threads to be used for processing. A value
#' below 1 indicates that all available threads will be utilized. Refer to the details section for
#' more information.
#' @param ... Additional arguments for writing raster outputs.
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
proximity <- function(x, samples_xy, radius_km = 200, num_threads = -1) {

    # check samples and histograms
    if (.is_mat(samples_xy)) {
        samples_xy <- .check_mat(samples_xy)
    } else {
        stop("'samples_xy' must be a matrix or an object convertibe to matrix.")
    }

    # check terra is available
    .check_pkgs("terra")
    # get the raster layers
    x <- .check_rast(x)

    # correction scale for long-lat CRS
    correction <- ifelse(.is_lonlat(x), 100000, 1)

    tryCatch(
        {
            output <- terra::app(
                x = x,
                fun = density_cpp,
                xy = samples_xy[, 1:2],
                radius_km = radius_km,
                scale = correction,
                nunum_threads = num_threads
                # na.rm = TRUE,
                # ...
            )
        },
        error = function(cond) {
            stop("Sample density calculation failed!\n", cond)
        }
    )

    return(
        output
    )
}

