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
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from \code{\link[terra]{predict}.
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
proximity <- function(x, samples_xy, radius_km = 200, num_threads = -1, ...) {

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
            output <- terra::interpolate(
                object = x,
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
    # check for NAs
    has_na <- anyNA(newdata)
    nr <- nrow(newdata)

    if (has_na) {
        idx <- which(stats::complete.cases(newdata))
        out <- rep(NaN, nr)
        # if all NA, return NaN vector
        if (!length(idx)) return(out)
        # subset the complete data
        dat <- as.matrix(newdata[idx, ])
    } else {
        dat <- as.matrix(newdata)
    }

    tryCatch(
        {
            pcount <- density_cpp(
                rast = dat,
                ...
            )
        },
        error = function(cond) {
            message("Error: the proximity C++ function faild, returning -0.02!")
            # return error values -0.02
            return(
                rep(-0.02, nr)
            )
        }
    )

    # sort out possible NAs
    if (has_na) {
        out[idx] <- pcount
        return(
            out
        )
    } else {
        return(
            pcount
        )
    }
}

