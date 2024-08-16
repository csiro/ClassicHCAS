#' The HCAS histogram
#'
#' The HCAS histogram calculation based on pair-point densities.
#'
#' You need to have OpenMP on your system to be able to benefit from parallel
#' processing to seed up the computations.
#'
#' @param observed The data.frame or matrix of observed remote sensing variables.
#' @param predicted The data.frame or matrix of predicted remote sensing variables.
#' @param samples Matrix or df. With the XY or XY and RS and ENV values.
#' @param within_km Numeric. Search radius in kilometers considering samples.
#' @param bin_width Numeric. The bin width of the histogram
#' @param bin_num Int. Number of bin for histogram. Keep it at the default 650..
#' @param num_threads Int. Number of CPU threads for processing...
#' @param filename Char (optional). The output file name for the .text file.
#'
#' @seealso \code{\link{benchmark}}
#'
#' @return matrix
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' obs_dat <- read.csv("")
#' prd_dat <- read.csv("")
#'
#'
#' }
histogram <- function(
        observed,
        predicted,
        samples,
        within_km = 1000,
        bin_width = 0.05,
        bin_num = 650,
        num_threads = parallel::detectCores() - 1,
        filename = "",
        source_code) {

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
        stop("There's NA in the extracted observed, predicted or xy samples!")

    if(any(colnames(predicted) != colnames(observed))) {
        warning("The names\\order of observed and predicted datasets doesn't match!\n")
        cat("Observed: ", colnames(observed), "\n")
        cat("Predicted:", colnames(predicted), "\n")
    }

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

    # write down the histogram
    if (nchar(filename) > 0) {
        filename <- ifelse(grepl(".txt", filename), filename, paste0(filename, ".txt"))
        tryCatch(
            {
                write.table(
                    out_table,
                    file = filename,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE
                )
            },
            error = function(cond) {
                stop("Writing histogram file failed!\n", cond)
            }
        )
    }

    return(
        out_table
    )
}

