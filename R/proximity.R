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
#' @param scale_factor Numeric. A scaling factor to correct distance calculations by converting degrees to
#' meters in unprojected coordinate systems. The default is 100,000 for an average conversion in Australia.
#' On the equator, this factor is approximately 111,235, and it varies with latitude according to a cosine
#' function.
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
proximity <- function(
        x,
        samples_xy,
        radius_km = 200,
        scale_factor = 100000,
        num_threads = -1,
        ...) {

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
    correction <- ifelse(.is_lonlat(x), scale_factor, 1)

    tryCatch(
        {
            output <- terra::interpolate(
                object = x[[1]],
                model = list(),
                fun = proxy_count,
                xy = samples_xy[, 1:2],
                radius_km = radius_km,
                scale = correction,
                num_threads = num_threads,
                ...
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


# wrapper function for density_cpp
proxy_count <- function(model, newdata, ...) {
    nr <- nrow(newdata)

    tryCatch(
        {
            pcount <- density_cpp(
                rast = as.matrix(newdata),
                ...
            )
        },
        error = function(cond) {
            message("Error: the proximity C++ function faild, returning -2!")
            # return error values -0.02
            return(
                rep(-2, nr)
            )
        }
    )

    return(pcount)
}

