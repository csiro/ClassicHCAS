#' Benchmarking target points
#'
#' The HCAS (Habitat Condition Assessment System) benchmarking function evaluates habitat
#' condition by comparing observed and predicted remote sensing (RS) variables. It integrates
#' multiple RS data layers to provide a comprehensive assessment of habitat quality and changes
#' over time. This function is designed to help researchers and conservationists quantify the
#' impacts of environmental changes and management interventions on habitat condition.
#'
#' Ensure that the order of RS variables is consistent between observed and predicted inputs
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
#' @param data A matrix, SpatRaster (from the \pkg{terra} package), or data.frame containing the input data.
#' The data must be organized in the following order: \strong{x}, \strong{y}, \strong{observed-RS},
#' \strong{predicted-RS} variables. If using a matrix or data.frame, ensure that the variables are in the
#' correct order. For a SpatRaster object, if it does not include x and y coordinates (or longitude and
#' latitude), you should set \code{add_xy = TRUE}. For more see the details section.
#' @param samples A matrix or data.frame containing x, y, observed-RS, and predicted-RS values
#' (in that specific order) for benchmark samples (also known as reference sites). If the
#' \code{data} argument is a SpatRaster, you can provide only the x and y coordinates of the
#' benchmark samples. In this case, the corresponding values will be extracted from the raster.
#' @param histogram A matrix or \strong{histo} object of normalised HCAS histogram
#' see \code{\link{histogram}} and \code{\link{normalise}}).
#' @param add_xy Logical. Add x and y to \code{data} (raster layers). If x and y are already
#' available as the first and second layers set \code{add_xy = FALSE}, otherwise result will not
#' be correct. For matrix \code{data} input this is ignored.
#' @param xy_stats A vector, mean and standard deviation of coordinates for centre and
#' scaling the coordinate to use as a penalty. The order should be: mean(x), mean(y), sd(x), sd(y)
#' @param xy_penalty Numeric. The spatial distance penalty value for selecting benchmark points.
#' The higher the value the more penalise the distant location will be. The value 0 means no penalty.
#' @param radius_km Numeric. Search radius in kilometers for considering benchmark samples.
#' @param k_pred Integer. Number of nearest ENV/predicted RS samples to take.
#' @param k_obs Integer. Number of nearest observed RS sample to takes.
#' @param bin_width Numeric. Specifies the bin width of the histogram. If \code{histogram} is
#' a \strong{histo} object, this value can be read from its attributes and may be left \code{NULL}.
#' The bin width must be consistent between the histogram creation and the benchmarking step to
#' ensure condition is accurately calculated.
#' @param interpolate Logical. Should interpolate the histogram?
#' @param offset Integer. Specifies the number of histogram bins that were ignored during
#' normalization (see \code{\link{clean}}).
#' @param confidence Numeric. The confidence value for LDC methods. See details below..
#' @param lambda Numeric. The lambda param for LDC Cauchy weighting...
#' @param make_su Logical. To make the uncertainty map or not.
#' @param num_threads Integer. Specifies the number of CPU threads to be used for processing. A value
#' below 1 indicates that all available threads will be utilized. Refer to the details section for
#' more information.
#' @param filename Char (optional). The output file name for raster outputs.
#' @param wopt list (optional). The output \code{\link[terra]{writeRaster}} options.
#' @param overwrite Logical. Whether to overwrite the output raster.
#'
#' @seealso \code{\link{histogram}}, \code{\link{normalise}}, and \code{\link{calibrate}}
#'
#' @return A matrix or SpatRaster
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' obs_dat <- read.csv("")
#' prd_dat <- read.csv("")
#' histo <- read.table("histogram.txt")
#'
#'
#' }
benchmark <- function(
        data,
        samples,
        histogram,
        add_xy = FALSE,
        xy_stats = c(0, 0, 1, 1),
        xy_penalty = 0.0,
        radius_km = 200,
        k_pred = 50,
        k_obs = 20,
        bin_width = NULL,
        interpolate = TRUE,
        offset = 0,
        confidence = 0.5,
        lambda = 2.0,
        exclude_slef = TRUE,
        make_su = FALSE,
        num_threads = -1,
        filename = "",
        wopt = list(datatype = "FLT4S", memfrac = 0.5),
        overwrite = TRUE) {

    # check samples and histograms
    if (.is_mat(samples)) {
        samples <- .check_mat(samples)
    } else {
        stop("'samples' must be a matrix or an object convertibe to matrix")
    }
    if (.is_mat(histogram)) {
        histogram <- .check_mat(histogram)
    } else {
        stop("'histogram' must be a matrix or an object convertibe to matrix")
    }

    if (nrow(histogram) != ncol(histogram)) {
        warning("Historgram dimensions are not the same!\n")
    }
    # cat("Histogram dimension:", dim(histogram), "\n")

    if (is(histogram, "histo")) {
        # check for histo offset consistency
        if (is.null(offset)) {
            offset <- attributes(histogram)$offset
        } else {
            if (offset != attributes(histogram)$offset) {
                warning("The supplied 'offset` is different from the arrtibute(histogram)$offset from the input.")
            }
        }
        # check for histo bin_width consistency
        if (is.null(bin_width)) {
            bin_width <- attributes(histogram)$bin.width
        } else {
            if (bin_width != attributes(histogram)$bin.width) {
                warning("The supplied 'bin_width` is different from the arrtibute(histogram)$bin.width from the input.")
            }
        }
    }

    # interpolate histogram
    if (interpolate) {
        histogram <- terra::as.matrix(
            x = terra::disagg(
                x = terra::rast(histogram),
                fact = 2,
                method = "bilinear"
            ),
            wide = TRUE
        )
    }
    # get the bin number after interpolation
    bin_num <- min(dim(histogram))

    # correction scale for long-lat CRS
    correction <- ifelse(.is_lonlat(data), 100000, 1)

    if (.is_mat(data)) {
        tryCatch(
            {
                # output <- bench_cpp(
                #     rast_stack = data,
                #     sample_vals = samples,
                #     histogram = histogram,
                #     xy_stats = xy_stats,
                #     xy_penalty = xy_penalty,
                #     within_km = radius_km,
                #     num_vars = 0, # keep it 0 for now to be calculated from data
                #     scale = correction,
                #     bin_width = ifelse(interpolate, bin_width / 2, bin_width),
                #     bin_num = bin_num,
                #     offset = ifelse(interpolate, offset * 2, offset),
                #     k_env = k_pred,
                #     k_rs = k_obs,
                #     confidence = confidence,
                #     lambda = lambda,
                #     exclude_slef = exclude_slef,
                #     make_su = make_su,
                #     num_threads = num_threads
                # )

                # change the class of samples to work with predict.hcas
                class(samples) <- c("hcas", "matrix", "array")

                output <- predict(
                    object = samples,
                    newdata = data,
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
                    make_su = make_su,
                    num_threads = num_threads
                )
            },
            error = function(cond) {
                stop("HCAS benchmarking failed!\n", cond)
            }
        )
    } else if (.is_rast(data)) {
        # get the raster layers
        data <- .check_rast(data)

        # add xy to the stack
        if (add_xy) {
            tryCatch(
                {
                    # NOTE: the order of files/layers matters here
                    # change the name to avoid error in in terra::predict
                    data <- c(
                        setNames(terra::init(predicted[[1]], fun = "x"), "x"),
                        setNames(terra::init(predicted[[1]], fun = "y"), "y"),
                        data
                    )
                },
                error = function(cond) {
                    stop("Failed to build xy for raster stack!\n", cond)
                }
            )
        }

        # sample extraction if needed
        if (ncol(samples) == 2) {
            cat("Extracting sample values...\n")
            samples <- as.matrix(
                terra::extract(data, samples, ID = FALSE)
            )
        } else if (ncol(samples) == terra::nlyr(data)) {
            cat("Samples with raster values are provided!\n")
        } else{
            # this should include rows/cols columns as well
            stop("Samples must either be coordinate values exclusively or include all raster layer values.")
        }

        # change the class of samples to work with terra::predict
        class(samples) <- c("hcas", "matrix", "array")

        tryCatch(
            {
                output <- terra::predict(
                    object = data,
                    model = samples,
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
                    make_su = make_su,
                    num_threads = num_threads,
                    na.rm = FALSE,
                    filename = filename,
                    overwrite = overwrite,
                    wopt = wopt
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
#' @export
#' @method predict hcas
predict.hcas <- function(object, newdata, make_su, ...){
    # check for NAs
    has_na <- anyNA(newdata)
    # number of output columns; TRUE/FALSE + 1
    nc <- make_su + 1
    nr <- nrow(newdata)

    if (has_na) {
        idx <- which(complete.cases(newdata))
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

    tryCatch(
        {
            # the HCAS C++ function
            hcas_cond <- bench_cpp(
                rast_stack = dat,
                sample_vals = as.matrix(object),
                make_su = make_su,
                ...
            )
        },
        error = function(cond) {
            message("Error: the benchmarking C++ function faild, returning -0.02!")
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

