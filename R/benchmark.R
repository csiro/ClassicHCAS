#' Benchmarking target points
#'
#' The HCAS benchmarking function to calculate habitat condition based on observed and
#' predicted remote sensing (RS) variables.
#'
#' Note that the order of RS variable must be the same for the observed and predicted
#' input data (in both raster and matrix inputs). The values of observed and
#' predicted RS variables must be centred and scaled (before prediction).
#'
#' You need to have OpenMP on your system to be able to benefit from parallel
#' processing to seed up the computations.
#'
#' The default parameters are designed for Australia and might nor work elsewhere.
#'
#'
#' @param data A matrix or SpatRaster (or data.frame) of the input data. The data must come
#' in the order of: x, y, observed-RS, predicted-RS variables. This could be a matrix or
#' data.frame object with the correct order of input variables or a SpatRaster object
#' from \pkg{terra} package. If the SpatRaster object does not include x and y coordinates (or
#' longitude and latitude), you must set \code{add_xy = TRUE}. See details for more
#' information.
#' @param samples A matrix or data.frame with the XY or alternatively, XY-observed-predicted values with
#' this order.
#' @param histogram A matrix of HCAS histogram (see \code{\link{histogram}}).
#' @param add_xy Logical. Add x and y to \code{data} (raster layers). If x and y are already
#' available as the first and second layers set \code{add_xy = FALSE}, otherwise result will not
#' be correct. For matrix \code{data} input this is ignored.
#' @param xy_stats A vector, mean and standard deviation of coordinates for centre and
#' scaling the coordinate to use as a penalty. The order should be: mean(x), mean(y), sd(x), sd(y)
#' @param xy_penalty Numeric. The penalty value the higher the value the more penalise the
#' distant location will be. The value 0 means no penalty
#' @param radius_km Numeric. Search radius in kilometers for considering benchmark samples.
#' @param k_pred Integer. Number of nearest ENV/predicted RS samples to take.
#' @param k_obs Integer. Number of nearest observed RS sample to takes.
#' @param bin_width Numeric. The bin width of the histogram
#' @param interpolate Logical. Should interpolate the histogram?
#' @param offset Int. Number of histogram bins to ignore...
#' @param confidence Numeric. The confidence value for LDC methods. See details below..
#' @param lambda Numeric. The lambda param for LDC Cauchy weighting...
#' @param make_su Logical. To make the uncertainty map or not.
#' @param num_threads Int. Number of CPU threads for processing...
#' @param filename Char (optional). The output file name for raster outputs
#' @param wopt list (optional). The output \code{\link[terra]{writeRaster}} options
#' @param overwrite Logical. Whether to overwrite the output raster.
#'
#' @seealso \code{\link{histogram}}
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
        bin_width = 0.05,
        interpolate = TRUE,
        offset = 0,
        confidence = 0.5,
        lambda = 2.0,
        exclude_slef = TRUE,
        make_su = FALSE,
        num_threads = parallel::detectCores() - 1,
        filename = "",
        wopt = list(datatype = "FLT4S", memfrac = 0.5),
        overwrite = TRUE) {

    # check samples and histograms
    if (.is_mat(samples)) {
        samples <- .check_mat(samples)
    } else {
        stop("'samples' must be a matrix or convertibe to matrix")
    }
    if (.is_mat(histogram)) {
        histogram <- .check_mat(histogram)
    } else {
        stop("'histogram' must be a matrix or convertibe to matrix")
    }

    if (nrow(histogram) != ncol(histogram)) {
        warning("Historgram dimensions are not the same!\n")
    }
    cat("Histogram dimension:", dim(histogram), "\n")
    # get the bin number after interpolation
    bin_num <- min(dim(histogram))
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

    # correction scale for long-lat CRS
    correction <- ifelse(.is_lonlat(data), 100000, 1)

    if (.is_mat(data)) {
        tryCatch(
            {
                output <- bench_cpp(
                    rast_stack = data,
                    sample_vals = samples,
                    histogram = histogram,
                    xy_stats = xy_stats,
                    xy_penalty = xy_penalty,
                    within_km = radius_km,
                    num_vars = num_vars,
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
                    num_vars = num_vars,
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

    }

    return(
        output
    )
}


# cat("Dimension of observed raster: ", dim(observed), "\n")
# cat("Dimension of predicted raster:", dim(predicted), "\n")
# len_non_na <- as.numeric(terra::global(predicted[[1]], fun = "notNA"))
# cat("Number of non-NA raster cells:", len_non_na, "\n")
# cat("The raster tile ")
# print(terra::ext(predicted[[1]]))
# if (terra::nlyr(observed) == terra::nlyr(predicted)) {
#     num_vars <- terra::nlyr(observed)
# } else {
#     stop("Number of observed and predicted raster layers are not equal!\n")
# }
#
# # check raster dimensions...
# if (any(dim(observed) != dim(predicted))) {
#     cat("WARNING: dimension of observed and predicted raster layers are different!\n")
# }
# if (any(c("x", "y", names(observed), names(predicted)) != colnames(samples))) {
#     cat("\nWARNING: Raster layers and samples column names do not match!\n")
#     cat("Observed layers: ", names(observed), "\n")
#     cat("Predicted layers:", names(predicted), "\n")
# }

# # if the raster (tile) is empty won't go ahead
# if (len_non_na < 1) {
#     if (make_su) {
#         # make two for SU; use first in case only one raster is used
#         out_map <- c(predicted[[1]], predicted[[1]])
#     } else {
#         out_map <- predicted[[1]]
#     }
#     names(out_map) <- paste0("lyr", seq_len(terra::nlyr(out_map)))
#     if (filename != "") {
#         terra::writeRaster(
#             x = out_map,
#             filename = filename,
#             overwrite = overwrite,
#             datatype = "FLT4S",
#             memfrac = 0.5
#         )
#     }
#
#     return(
#         out_map
#     )
# }


# the generic HCAS prediction function for C++ integration with terra
# NOTE: the na.rm arg in terra::predict doesn't provide the correct mask
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

