#' Calibrate habitat condition output
#'
#' This function calibrates the HCAS habitat condition values and scales them between 0 and 1
#' using a Spline function.
#'
#' @param x A SpatRaster, matrix, data.frame, or vector containing HCAS habitat condition values from
#' the benchmarking function (see \code{\link{benchmark}}). If \code{x} is a matrix or data.frame,
#' the calibration will be applied to all columns.
#' @param x_values Numeric vector of uncalibrated condition values. It is recommended that \code{x_values}
#' cover the full range of values, including the minimum and maximum of the raw condition data.
#' @param y_values Numeric vector of calibrated target condition value corresponding to \code{x_values}.
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.
#'
#' @seealso \code{\link{benchmark}}
#'
#' @return A matrix, SpatRaster or vector, depending on the inputs.
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#'
#'
#' }
calibrate <- function(x, x_values, y_values = c(0, 0.101, 0.944, 1), ...) {
    # some initial checks
    if (!methods::is(x_values, "numeric")) {
        stop("The 'x_values' must be a numeric vector.")
    }
    if (!methods::is(y_values, "numeric")) {
        stop("The 'y_values' must be a numeric vector.")
    }
    # the length of the x and y should be the same
    if (length(x_values) != length(y_values)) {
        stop("The length of 'x_values' and 'y_values' must be the same.")
    }

    # calibration
    if (.is_rast(x)) {
        # check terra is available
        .check_pkgs("terra")
        # check x
        x <- .check_rast(x)
        # calibrate condition raster
        out <- terra::app(
            x = x,
            fun = .calib,
            x_vals = x_values,
            y_vals = y_values,
            ...
        )
    } else {
        # if matrix or data.frame apply on all columns
        if (.is_mat(x)) {
            # check and convert to matrix
            x <- .check_mat(x)

            out <- apply(
                X = x,
                MARGIN = 2,
                FUN = .calib,
                x_vals = x_values,
                y_vals = y_values
            )
        } else {
            out <- .calib(
                x = x,
                x_vals = x_values,
                y_vals = y_values
            )
        }
    }

    return(out)
}


# a general function to calibrate condition
# NOTE: although this will repeat spline fitting when raster files are read in steps,
#   the spline function is deterministic and very fast, 0.005 seconds
.calib <- function(x, x_vals, y_vals) {
    # fitting smoothing spline
    spline_dat <- as.data.frame(
        stats::spline(
            x= x_vals,
            y = y_vals,
            method = "natural",
            xout = seq(0, x_vals[length(x_vals)], length.out = 100)
        )
    )

    return(
        spline_rescale(
            x = x,
            spline_x = spline_dat$x,
            spline_y = spline_dat$y
        )
    )
}

