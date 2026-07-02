#' Calibrate HCAS condition values
#'
#' Transforms raw, unscaled HCAS benchmarking values to a standard 0-1 habitat
#' condition scale using a monotonic spline.
#'
#' @details
#' \code{\link{benchmark}} returns an unscaled relative condition value. The
#' range and distribution of that value can differ among regions, data sources,
#' remote-sensing summaries, and model settings. \code{calibrate()} maps those
#' raw values to an interpretable 0-1 scale, where 0 represents a completely
#' degraded or removed state and 1 represents reference or near-natural
#' condition.
#'
#' Calibration uses paired \code{x_values} and \code{y_values}. The
#' \code{x_values} are raw condition values from the benchmark output, often
#' including the minimum, the median of highly modified sites, the median of
#' reference sites, and the maximum. The \code{y_values} are the desired
#' calibrated condition values for those knots, supplied from empirical evidence,
#' expert judgement, or another accepted condition scale. A monotonic spline is
#' fitted through the pairs, preserving the rank order of raw condition while
#' enforcing a smooth increasing transformation. Values outside the 0-1 range
#' after transformation are clipped to 0 or 1.
#'
#' \code{x_values} and \code{y_values} must have the same length. For a stable
#' calibration curve, \code{x_values} should be sorted from low to high and
#' should span the raw condition values in \code{x}. The default
#' \code{y_values} assumes four calibration knots and should be replaced when a
#' different number of \code{x_values} is supplied.
#'
#' @param x A \pkg{terra} \code{SpatRaster}, matrix, data.frame, or vector
#' containing raw HCAS condition values returned by \code{\link{benchmark}}. For
#' matrix and data.frame inputs, calibration is applied column-wise.
#' @param x_values Numeric vector of raw, uncalibrated condition values used as
#' calibration knots.
#' @param y_values Numeric vector of calibrated target values corresponding to
#' \code{x_values}. Values should usually be between 0 and 1.
#' @param ... Additional arguments passed to \code{\link[terra]{app}} when
#' calibrating raster outputs, such as \code{filename}, \code{overwrite}, or
#' \code{wopt}.
#'
#' @seealso \code{\link{benchmark}}
#'
#' @return A vector, matrix, or \pkg{terra} \code{SpatRaster}, depending on the
#' input.
#' @export
#'
#' @examples
#' library(ClassicHCAS)
#'
#' raw_condition <- c(0.000, 0.005, 0.020, 0.035, 0.050)
#'
#' calibrated <- calibrate(
#'     raw_condition,
#'     x_values = c(0.000, 0.005, 0.035, 0.050),
#'     y_values = c(0.0, 0.1, 0.9, 1.0)
#' )
#'
#' calibrated
calibrate <- function(
        x,
        x_values,
        y_values = c(0, 0.101, 0.944, 1),
        ...) {

    dots <- list(...)
    if ("interpolate" %in% names(dots)) {
        stop("The 'interpolate' argument is no longer supported by 'calibrate()'.")
    }

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
.calib <- function(x, x_vals, y_vals) {
    # the monotonic spline function
    f <- stats::splinefun(
        x = x_vals,
        y = y_vals,
        method = "monoH.FC"
    )

    out <- rep(NA_real_, length(x))
    ok <- is.finite(x)
    if (any(ok)) out[ok] <- f(x[ok])
    out[out > 1] <- 1
    out[out < 0] <- 0

    return(out)
}

