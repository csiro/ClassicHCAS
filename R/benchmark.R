#' Condition benchmarking of target points
#'
#' The HCAS (Habitat Condition Assessment System) benchmarking function evaluates habitat
#' condition by comparing observed and predicted remote sensing (RS) variables. It integrates
#' multiple RS data layers to provide a comprehensive assessment of habitat quality and changes
#' over time. This function is designed to help researchers and conservationists quantify the
#' impacts of environmental changes and management interventions on habitat condition.
#'
#' Ensure that the order of remote sensing variables is consistent between observed and predicted inputs
#' (for both raster and matrix formats). The RS variable values must be centered and scaled
#' prior to prediction. Failure to do so may result in variables with larger ranges having
#' disproportionate influence in the multi-dimensional distance calculations.
#'
#' Ensure that \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system
#' to take advantage of parallel processing and accelerate computations. While most systems
#' include OpenMP by default, you may need to load the appropriate module if you're using an HPC
#' system.
#'
#' Note that the default parameters are tailored for Australia and may not be suitable for other
#' regions.
#'
#' @inheritParams histogram
#' @param samples A matrix or data.frame containing x, y, observed-RS, and predicted-RS values
#' (in that specific order) for benchmark samples (also known as reference sites). If the
#' \code{data} argument is a SpatRaster, you can provide only the x and y coordinates of the
#' benchmark samples. In this case, the corresponding values will be extracted from the raster
#' layers. Consider extra time for sample value extraction in this case.
#' @param histogram A matrix or \strong{histo} object of normalised HCAS histogram
#' see \code{\link{histogram}} and \code{\link{normalise}}).
#' @param xy_stats A vector, mean and standard deviation of coordinates for centre and
#' scaling the coordinate to use as a penalty. The order should be: mean(x), mean(y), sd(x), sd(y)
#' @param xy_penalty Numeric. The spatial distance penalty value for selecting benchmark points.
#' The higher the value the more penalise the distant location will be. The value 0 means no penalty.
#' @param radius_km Numeric. Search radius in kilometers for considering benchmark samples.
#' @param k_pred Integer. Number of nearest predicted RS samples to take.
#' @param k_obs Integer. Number of nearest observed RS sample to takes.
#' @param bin_width Numeric. Specifies the bin width of the histogram. If \code{histogram} is
#' a \strong{histo} object, this value can be read from its attributes and may be left \code{NULL}.
#' The bin width must be consistent between the histogram creation and the benchmarking step to
#' ensure condition is accurately calculated.
#' @param interpolate Logical. Whether to interpolate the histogram for a smoother result.
#' @param offset Integer. Specifies the number of histogram bins that were ignored during
#' normalisation (see \code{\link{normalise}}). If \code{histogram} is a \strong{histo} object,
#' this value can be read from its attributes and may be set \code{NULL}. Similar to bin-width,
#' the \code{offset} must be consistent between the histogram normalisation and the benchmarking
#' step to ensure condition is accurately calculated.
#' @param confidence Numeric. The confidence value for LDC methods. See details below..
#' @param lambda Numeric. The lambda param for LDC Cauchy weighting...
#' @param exclude_slef Logical. To exclude a benchmark point from assessing itself.
#' @param drop_features Integer vector. Completely remove the RS variable from the benchmarking process. For
#' consistency, it is recommended to exclude the same variables used in the histogram step; unless
#' you have a specific reason not to.
#' @param make_su Logical. To make the uncertainty map or not.
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.
#'
#' @seealso \code{\link{histogram}}, \code{\link{normalise}}, and \code{\link{calibrate}}
#'
#' @return A matrix or SpatRaster, depending on the inputs.
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#'
#'
#' }
benchmark <- function(
        data,
        samples,
        histogram,
        xy_stats = c(0, 0, 1, 1),
        xy_penalty = 0.0,
        radius_km = 200,
        scale_factor = 100000,
        k_pred = 50,
        k_obs = 20,
        bin_width = NULL,
        interpolate = TRUE,
        offset = 0,
        confidence = 0.5,
        lambda = 2.0,
        exclude_slef = TRUE,
        drop_features = NULL,
        make_su = FALSE,
        num_threads = -1,
        ...) {

    # check k_pred and k_obs
    if (k_pred < k_obs) stop("'k_obs' must be smaller or equal to 'k_pred'.")
    # check samples and histograms
    samples <- if (.is_mat(samples)) .check_mat(samples) else stop("'samples' must be a matrix or convertible to one.")
    histogram <- if (.is_mat(histogram)) .check_mat(histogram) else stop("'histogram' must be a matrix or convertible to one.")
    if (nrow(histogram) != ncol(histogram)) warning("Historgram dimensions are not the same!\n")

    if (methods::is(histogram, "histo")) {
        # check for histo bin_width consistency
        if (is.null(bin_width)) {
            bin_width <- attributes(histogram)$bin.width
        } else {
            if (bin_width != attributes(histogram)$bin.width) {
                warning("Provided 'bin_width' differs from histogram attribute.")
            }
        }
        # check for histo offset consistency
        if (is.null(offset)) {
            offset <- attributes(histogram)$offset
        } else {
            if (offset != attributes(histogram)$offset) {
                warning("Provided 'offset' differs from histogram attribute.")
            }
        }
    }

    # interpolate histogram
    if (interpolate) {
        histogram <- terra::as.matrix(
            terra::disagg(terra::rast(histogram), fact = 2, method = "bilinear"),
            wide = TRUE
        )
    }
    # get the bin number after interpolation
    bin_num <- min(dim(histogram))

    exclude_var <- NULL
    # drop features from calculation if requested
    if (length(drop_features)) {
        n_vars <- (ncol(samples) - 2) / 2
        exclude_var <- c(drop_features + 2, drop_features + 2 + n_vars)
    }

    if (.is_mat(data)) {
        # check and convert to matrix
        data <- .check_mat(data)
        # correction scale for long-lat CRS
        correction <- ifelse(.is_lonlat(data), scale_factor, 1)

        if (ncol(samples) != ncol(data)) {
            stop("Samples must include all raster values (matching column count with 'data').")
        }

        # remove the features from the reference samples as well
        if (!is.null(exclude_var)) samples[, exclude_var] <- 0

        tryCatch(
            {
                output <- benchmarking(
                    model = list(),
                    newdata = data, # rast_stack arg
                    sample_vals = samples,
                    histogram = histogram,
                    xy_stats = xy_stats,
                    xy_penalty = xy_penalty,
                    within_km = radius_km,
                    num_vars = 0, # keep it 0 for now to be calculated from data
                    scale = correction,
                    bin_width = ifelse(interpolate, bin_width / 2, bin_width),
                    bin_num = bin_num,
                    offset = ifelse(interpolate, offset * 2, offset),
                    k_env = k_pred,
                    k_rs = k_obs,
                    confidence = confidence,
                    lambda = lambda,
                    exclude_slef = exclude_slef,
                    drop = exclude_var,
                    make_su = make_su,
                    num_threads = num_threads
                )
            },
            error = function(cond) {
                stop("HCAS benchmarking failed!\n", cond)
            }
        )
    } else if (.is_rast(data)) {
        # check terra is available
        .check_pkgs("terra")
        # check and convert to SpatRaster object
        data <- .check_rast(data)

        # correction scale for long-lat CRS
        correction <- ifelse(.is_lonlat(data), scale_factor, 1)

        # sample extraction if needed
        if (ncol(samples) == 2) {
            cat("Extracting sample values...\n")
            samples <- cbind(samples, as.matrix(terra::extract(data, samples)))
        } else if ((ncol(samples) - 2) != terra::nlyr(data)) {
            stop("Sample feature count does not match number of raster layers.")
        }

        # remove the features from the reference samples as well
        if (length(drop_features)) {
            samples[, exclude_var] <- 0
        }

        tryCatch(
            {
                output <- terra::interpolate(
                    object = data,
                    model = list(),
                    fun = benchmarking,
                    sample_vals = samples,
                    histogram = histogram,
                    xy_stats = xy_stats,
                    xy_penalty = xy_penalty,
                    within_km = radius_km,
                    num_vars = 0, # keep it 0 for now to be calculated from data
                    scale = correction,
                    bin_width = ifelse(interpolate, bin_width / 2, bin_width),
                    bin_num = bin_num,
                    offset = ifelse(interpolate, offset * 2, offset),
                    k_env = k_pred,
                    k_rs = k_obs,
                    confidence = confidence,
                    lambda = lambda,
                    exclude_slef = exclude_slef,
                    drop = exclude_var,
                    make_su = make_su,
                    num_threads = num_threads,
                    ...
                )
            },
            error = function(cond) {
                stop("HCAS benchmarking failed!\n", cond)
            }
        )

    } else {
        stop("The 'data' must be raster or a matrix, or convertiable object to these classes.")
    }

    return(
        output
    )
}


# the generic HCAS prediction function for C++ integration with terra
# NOTE: the na.rm arg in terra::predict doesn't provide the correct mask
#   and keep model an empty list
benchmarking <- function(model, newdata, make_su, drop = NULL, ...){
    # check for NAs
    has_na <- anyNA(newdata)
    # number of output columns; TRUE/FALSE + 1
    nc <- make_su + 1
    nr <- nrow(newdata)

    if (has_na) {
        idx <- which(stats::complete.cases(newdata))
        out <- matrix(NaN, nrow = nr, ncol = nc)
        # name the raster layers
        colnames(out) <- c("condition", "su")[1:nc]
        # if all NA, return NaN vector
        if (!length(idx)) return(out)
        # subset the complete data
        dat <- as.matrix(newdata[idx, ])
    } else {
        dat <- as.matrix(newdata)
    }

    # if drop_features are provided exclude them from features
    if (length(drop)) dat[, drop] <- 0

    tryCatch(
        {
            # the HCAS C++ function
            hcas_cond <- bench_cpp(
                rast_stack = dat,
                make_su = make_su,
                ...
            )
        },
        error = function(cond) {
            message("Benchmarking C++ function failed. Returning -0.02 for all cells.")
            # return error values -0.02
            return(
                matrix(-0.02, nrow = nr, ncol = nc)
            )
        }
    )

    if (has_na) {
        out[idx, ] <- hcas_cond
        return(out)
    } else {
        # name the raster layers
        colnames(hcas_cond) <- c("condition", "su")[1:nc]
        return(hcas_cond)
    }
}

