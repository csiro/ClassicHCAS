#' Benchmarking
#'
#' The HCAS benchmarking function to calculate habitat condition based on observed and
#' predicted remote sensing variables.
#'
#' @param observed The data.frame or matrix of observed remote sensing variables.
#' @param predicted
#' @param samples
#' @param histogram
#' @param xy_stats
#' @param xy_penalty
#' @param within_km
#' @param k_env
#' @param k_rs
#' @param bin_width
#' @param interpolate
#' @param offset
#' @param confidence
#' @param lambda
#' @param make_su
#' @param force_compile
#' @param num_threads
#' @param filename
#' @param overwrite
#' @param source_code
#'
#' @return data.frame or
#' @export
#'
#' @examples
benchmarking <- function(
        observed,
        predicted,
        samples,
        histogram,
        xy_stats = c(0, 0, 1, 1),
        xy_penalty = 0.0,
        within_km = 200,
        k_env = 50,
        k_rs = 20,
        bin_width = 0.03,
        interpolate = TRUE,
        offset = 0,
        confidence = 0.5,
        lambda = 2.0,
        make_su = FALSE,
        force_compile = FALSE,
        num_threads = parallel::detectCores() - 1,
        filename = "",
        overwrite = TRUE,
        source_code) {

    # compile the C++ code
    # cat("Compiling C++ code...\n")
    # tryCatch(
    #     {
    #         Rcpp::sourceCpp(source_code, rebuild = force_compile)
    #     },
    #     error = function(cond) {
    #         stop("Compiling C++ file failed!\n", cond)
    #     }
    # )

    # get the raster layers
    observed <- get_rast(observed)
    predicted <- get_rast(predicted)

    # correction scale for long-lat CRS
    correction <- ifelse(terra::is.lonlat(predicted, perhaps = TRUE), 100000, 1)

    if (terra::nlyr(observed) == terra::nlyr(predicted)) {
        num_vars <- terra::nlyr(observed)
    } else {
        stop("Number of observed and predicted raster layers are not equal!\n")
    }

    if (nrow(histogram) != ncol(histogram)) {
        cat("WARNING: Historgram dimensions are not the same!\n")
    }
    cat("Histogram dimension:", dim(histogram), "\n")

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


    # check raster dimensions...
    if (any(dim(observed) != dim(predicted))) {
        cat("WARNING: dimension of observed and predicted raster layers are different!\n")
    }
    cat("Dimension of observed raster: ", dim(observed), "\n")
    cat("Dimension of predicted raster:", dim(predicted), "\n")
    len_non_na <- as.numeric(terra::global(predicted[[1]], fun = "notNA"))
    cat("Number of non-NA raster cells:", len_non_na, "\n")
    cat("The raster tile ")
    print(terra::ext(predicted[[1]]))

    # if the raster (tile) is empty won't go ahead
    if (len_non_na < 1) {
        if (make_su) {
            # make two for SU; use first in case only one raster is used
            out_map <- c(predicted[[1]], predicted[[1]])
        } else {
            out_map <- predicted[[1]]
        }
        names(out_map) <- paste0("lyr", seq_len(terra::nlyr(out_map)))
        if (filename != "") {
            terra::writeRaster(
                x = out_map,
                filename = filename,
                overwrite = overwrite,
                datatype = "FLT4S",
                memfrac = 0.5
            )
        }

        return(
            out_map
        )
    }

    if (any(c("x", "y", names(observed), names(predicted)) != colnames(samples))) {
        cat("\nWARNING: Raster layers and samples column names do not match!\n")
        cat("Observed layers: ", names(observed), "\n")
        cat("Predicted layers:", names(predicted), "\n")
    }
    # make the stack
    cat("\nBuilding raster stack...\n")
    tryCatch(
        {
            # NOTE: the order of files/layers matters here
            # change the name to avoid error in in terra::predict
            rast_stack <- c(
                setNames(terra::init(predicted[[1]], fun = "x"), "x"),
                setNames(terra::init(predicted[[1]], fun = "y"), "y"),
                setNames(observed, paste0("RS", 1:num_vars)),
                setNames(predicted, paste0("EN", 1:num_vars))
            )
        },
        error = function(cond) {
            stop("Failed to build raster stack!\n", cond)
        }
    )


    # sample extraction if needed
    if (ncol(samples) == 2) {
        cat("Extracting sample values...\n")
        samples <- as.matrix(
            terra::extract(rast_stack, samples, ID = FALSE)
        )
    } else if (ncol(samples) == 2*num_vars+2) {
        cat("Samples with raster values are provided!\n")
        if (!is(samples, "matrix"))
            samples <- as.matrix(samples)
    } else{
        # this should include rows/cols columns as well
        stop("Samples must either be coordinate values exclusively or include all raster layer values.")
    }

    cat("Reference samples dimension:", dim(samples), "\n")

    # change the class of samples
    class(samples) <- c("hcas", "matrix", "array")

    cat("\nComputing HCAS benchmarking...\n")
    tryCatch(
        {
            out_map <- terra::predict(
                object = rast_stack,
                model = samples,
                histogram = histogram,
                xy_stats = xy_stats,
                xy_penalty = xy_penalty,
                within_km = within_km,
                num_vars = num_vars,
                scale = correction,
                bin_width = ifelse(interpolate, bin_width / 2, bin_width),
                bin_num = bin_num,
                offset = ifelse(interpolate, offset * 2, offset),
                k_env = k_env,
                k_rs = k_rs,
                confidence = confidence,
                lambda = lambda,
                num_threads = num_threads,
                make_su = make_su,
                na.rm = FALSE,
                filename = filename,
                overwrite = overwrite,
                wopt = list(datatype = "FLT4S", memfrac = 0.5)
            )
        },
        error = function(cond) {
            stop("HCAS benchmarking failed!\n", cond)
        }
    )
    cat("HCAS benchmarking was successful!\n")

    return(
        out_map
    )
}


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
            hcas_cond <- hcas_cpp(
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

