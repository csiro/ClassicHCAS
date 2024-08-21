#' The HCAS histogram
#'
#' The HCAS histogram calculation is based on pair-point densities. It is an integral part of HCAS,
#' designed to learn the expected observed remote sensing (RS) values from the
#' predicted RS values across a wide range of reference sites. Note that the reference
#' sites used for the \code{histogram} do not need to be the same as those used in
#' the \code{\link{benchmark}} function. See more in details.
#'
#' Ensure that the order of RS variables is consistent between observed and predicted inputs.
#' The RS variable values must be centered and scaled prior to prediction. Failure to do
#' so may result in variables with larger ranges having disproportionate influence in the
#' multi-dimensional distance calculations.
#'
#' Ensure that \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system
#' to take advantage of parallel processing and accelerate computations. While most systems
#' include OpenMP by default, you may need to load the appropriate module if you're using an HPC
#' system.
#'
#' @param observed A matrix or data.frame of observed remote sensing variables.
#' @param predicted A matrix or data.frame of predicted remote sensing variables.
#' Notice that the order or predicted RS variables must the same as observed ones.
#' @param samples_xy A matrix or data.frame of x and y (longitude and latitude) of the
#' reference points used for observed and predicted RS variiables.
#' @param radius_km Numeric. Search radius in kilometers for considering reference samples
#' in creating histogram.
#' @param bin_width Numeric. Specifies the bin width of the histogram. Finding an optimal bin width
#' may require some experimentation to achieve the best results.
#' @param bin_num Integer. Specifies the number of bins for the histogram. It is generally recommended
#' to use the default value of 650. Adjusting the \code{bin_width} is often more effective than changing
#' \code{bin_num}.
#' @param num_threads Int. Number of CPU threads for processing. See details.
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
        samples_xy,
        radius_km = 1000,
        bin_width = 0.05,
        bin_num = 650,
        num_threads = parallel::detectCores() - 1,
        filename = "",
        source_code) {

    # check samples_xy
    if (.is_mat(samples_xy)) {
        samples_xy <- .check_mat(samples_xy)
    } else {
        stop("'samples_xy' must be a matrix or an object convertibe to matrix")
    }

    if (ncol(samples_xy) != 2)
        stop("'samples_xy' must be a data.frame or matrix with exactly two columns: one for longitude (x) and one for latitude (y).")

    # correction scale for long-lat CRS
    correction <- ifelse(.is_lonlat(samples_xy), 100000, 1)

    # check for predicted variables
    if (.is_mat(predicted)) {
        predicted <- .check_mat(predicted)
    } else if (.is_rast(predicted)) {
        # extract values
        predicted <- as.matrix(
            terra::extract(predicted, samples_xy, ID = FALSE)
        )
    }  else {
        stop("The 'predicted' must be raster or a matrix, or convertiable object to these classes.")
    }

    # check for observed variables
    if (.is_mat(observed)) {
        observed <- .check_mat(observed)
    } else if (.is_rast(observed)) {
        # extract values
        observed <- as.matrix(
            terra::extract(observed, samples_xy, ID = FALSE)
        )
    }  else {
        stop("The 'observed' must be raster or a matrix, or convertiable object to these classes.")
    }

    cat("\nRS data dimensions: ", dim(observed), "\n")
    cat("ENV data dimensions:", dim(predicted), "\n")

    # some error checking
    if(any(dim(predicted) != dim(observed)))
        stop("Dimensions of RS and ENV datasets doesn't match!")

    if(nrow(predicted) != nrow(samples_xy))
        stop("Number of rows of rasters values and reference samples_xy doesn't match!")

    if(any(anyNA(predicted), anyNA(observed), anyNA(samples_xy)))
        stop("There's NA in the extracted observed, predicted or xy samples_xy!")

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
                samples_xy = samples_xy,
                within_km = radius_km,
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

