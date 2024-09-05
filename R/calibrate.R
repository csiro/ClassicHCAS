#' Calibrate habitat condition output
#'
#' This function calibrates the HCAS habitat condition and scales the values between
#' 0 and 1.
#'
#' @param x A SpatRaster, matrix, data.frame, or vector containing HCAS habitat condition values from
#' the benchmarking function (see \code{\link{benchmark}}). If \code{x} is a matrix or data.frame,
#' the calibration will be applied to all columns.
#' @param low Numeric. The median uncalibrated condition value for highly modified areas.
#' @param mid Numeric. The median uncalibrated condition value for natural reference areas.
#' @param high Numeric. The median uncalibrated condition value for \strong{benchmark} reference
#' sites used in the \code{\link{benchmark}} function.
#' @param max Numeric. The median uncalibrated condition value for the entire area of interest, based
#' on the maximum values in the raster or matrix.
#' @param target1 Numeric. The calibrated target value for highly modified areas. The default value of
#' 0.1 is based on Australian data.
#' @param target2 Numeric. The calibrated target value for natural reference areas. The default value
#' of 0.944 is based on Australian data.
#' @param method Character. Specifies the calibration method, either "linear" or "spline". The spline
#' method produces a smoother map and histogram of calibrated condition values. The default is "spline".
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from \code{\link[terra]{predict}.
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
calibrate <- function(x, low, mid, high, max,
                      target1 = 0.1, target2 = 0.944,
                      method = "spline",
                      ...) {

    method <- match.arg(method, c("spline", "linear"))
    linear <- ifelse(method == "linear", TRUE, FALSE)

    if (.is_rast(x)) {
        # check terra is available
        .check_pkgs("terra")
        # check x
        x <- .check_rast(x)
        # calibrate condition raster
        out <- terra::app(
            x = x,
            fun = .calib,
            l = low,
            m = mid,
            h = high,
            mx = max,
            t1 = target1,
            t2 = target2,
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
                l = low,
                m = mid,
                h = high,
                mx = max,
                t1 = target1,
                t2 = target2,
                linear = linear
            )
        } else {
            out <- .calib(
                x = x,
                l = low,
                m = mid,
                h = high,
                mx = max,
                t1 = target1,
                t2 = target2,
                linear = linear
            )
        }
    }

    return(out)
}


# a general function to calibrate condition
# NOTE: although this will repeat spline fitting when raster files are read in steps,
# the spline function is deterministic and very fast, 0.005 seconds
.calib <- function(x, l, m, h, mx, t1, t2, linear = FALSE) {
    # linear or spline
    if (linear) {
        return(
            linear_rescale(
                x = x,
                low_point = l,
                high_point = m,
                max_point = h,
                low_target = t1,
                high_target = t2
            )
        )
    } else {
        # spline data
        dat <- data.frame(
            x = c(0, l, m, h, mx),
            y = c(0, t1, t2, 1, 1)
        )
        # fitting smoothing spline
        spline_dat <- as.data.frame(
            stats::spline(
                x= dat$x,
                y = dat$y,
                method = "natural",
                xout = seq(0, mx, length.out = 100)
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
}

