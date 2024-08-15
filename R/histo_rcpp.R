#' Histogram
#'
#' The HCAS histogram calcualtion based on pair-point densities.
#'
#' @param observed
#' @param predicted
#' @param samples
#' @param within_km
#' @param bin_width
#' @param bin_num
#' @param force_compile
#' @param num_threads
#' @param output_file
#' @param source_code
#'
#' @return
#' @export
#'
#' @examples
histo_hcas <- function(
        observed,
        predicted,
        samples,
        within_km = 1000,
        bin_width = 0.05,
        bin_num = 650,
        force_compile = FALSE,
        num_threads = parallel::detectCores() - 1,
        output_file,
        source_code) {

    # compile the C++ code
    cat("Compiling C++ code...\n")
    tryCatch(
        {
            Rcpp::sourceCpp(source_code, rebuild = force_compile)
        },
        error = function(cond) {
            stop("Compiling C++ file failed!\n")
        }
    )

    # check samples
    if (is(samples, "data.frame")) {
        samples <- as.matrix(samples)
    } else if (!is(samples, "matrix")) {
        stop("Samples must be a data.frame or matrix of the coordinate values.")
    }

    # correction scale for long-lat CRS
    correction <- ifelse(.is_lonlat(samples), 100000, 1)

    # check for predicted variables
    if(is(predicted, "data.frame")) {
        predicted <- as.matrix(predicted)
    } else if(is(predicted, "SpatRaster")) {
        # extract values
        predicted <- as.matrix(
            terra::extract(predicted, samples, ID = FALSE)
        )
    } else if(!is(predicted, "matrix")) {
        stop("'predicted' must be raster layers or a matrix of extracted ENV values.")
    }

    # check for observed variables
    if(is(observed, "data.frame")) {
        observed <- as.matrix(observed)
    } else if(is(observed, "SpatRaster")) {
        # extract values
        observed <- as.matrix(
            terra::extract(observed, samples, ID = FALSE)
        )
    } else if(!is(observed, "matrix")) {
        stop("'observed' must be raster layers or a matrix of extracted RS values.")
    }


    cat("\nRS data dimensions: ", dim(observed), "\n")
    cat("ENV data dimensions:", dim(predicted), "\n")

    # some error checking
    if(any(dim(predicted) != dim(observed)))
        stop("Dimensions of RS and ENV datasets doesn't match!")

    if(nrow(predicted) != nrow(samples))
        stop("Number of rows of rasters values and reference samples doesn't match!")

    if(any(anyNA(predicted), anyNA(observed), anyNA(samples)))
        stop("There's NA in the extracted RS\\ENV or samples!")

    if(any(colnames(predicted) != colnames(observed))) {
        cat("\nWARNING: The names\\order of OBS and MOD datasets doesn't match!\n")
        cat("Observed:", colnames(observed), "\n")
        cat("Modelled:", colnames(predicted), "\n")
    }

    cat("\nCalculating histogram...\n")
    # run histo calculate in C++
    tryCatch(
        {
            out_table <- histo_cpp(
                rs_vals = observed,
                pr_vals = predicted,
                samples_xy = samples,
                within_km = within_km,
                scale = correction,
                bin_width = bin_width,
                bin_num = bin_num,
                num_threads = num_threads
            )
        },
        error = function(cond) {
            stop("Histogram calculation failed!\n", cond)
        }
    )
    cat("Histogram is calculated successfully!\n")

    # write down the histogram
    tryCatch(
        {
            write.table(
                out_table,
                file = output_file,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
            )
        },
        error = function(cond) {
            stop("Writing histogram file failed!\n", cond)
        }
    )
    cat("Histogram is written successfully!\n")

    return(
        out_table
    )
}


# are sample points lat-log?
.is_lonlat <- function(x) {
    terra::is.lonlat(
        terra::vect(x),
        perhaps = TRUE,
        warn = FALSE
    )
}

