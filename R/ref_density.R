#' Build an HCAS reference density surface
#'
#' Builds the raw reference density, or probability, surface used by HCAS to
#' describe natural variation between predicted and observed remote-sensing (RS)
#' variables at reference sites.
#'
#' @details
#' In HCAS, predicted RS values represent the expected reference-condition signal
#' for a site, usually estimated from environmental covariates using
#' high-integrity reference ecosystems. Observed RS values are the actual Earth
#' observation summaries for the same sites. \code{ref_density()} compares
#' reference samples with one another and records how predicted-distance and
#' observed-distance co-vary under reference condition.
#'
#' For each pair of reference samples within \code{radius_km}, the function
#' calculates L1 (Manhattan) distances across all retained RS variables:
#'
#' \deqn{d(x, y) = \sum_{j = 1}^{m} |x_j - y_j|}
#'
#' where \eqn{m} is the number of RS variables. The predicted-distance and
#' observed-distance values are then added to a two-dimensional surface with
#' resolution controlled by \code{bin_width} and \code{bin_num}. The resulting
#' surface is a compact empirical summary of how much observed RS variation is
#' expected for a given predicted RS distance among intact or near-intact
#' reference samples.
#'
#' The output of this function is raw and should normally be passed to
#' \code{\link{normalise}} before use in \code{\link{benchmark}}. The
#' \code{bin.width} attribute is stored on the returned object so that the same
#' bin width can be reused during normalisation and benchmarking.
#'
#' Matrix and data.frame inputs must be ordered as \code{x}, \code{y}, predicted
#' RS variables, then observed RS variables. Raster inputs should contain the
#' predicted RS layers followed by the observed RS layers in the same variable
#' order; the reference coordinates are supplied separately through
#' \code{samples}. Predicted and observed variables should be centred and scaled
#' consistently before running this function so that no variable dominates the
#' multidimensional distance calculation because of its units or numeric range.
#'
#' In geographic coordinates, radius searches use a fast integer approximation.
#' Coordinates are stored in micro-degrees (\code{degree * 1000000}) and distance
#' is approximated by:
#'
#' \deqn{distance^2 \approx dlat^2 + (dlon \times \cos(lat_1))^2}
#'
#' where \eqn{\cos(lat_1)} is derived from the query latitude. This is efficient
#' for large analyses but introduces distortion over large areas. For high
#' accuracy at broad regional or continental radii, use a projected coordinate
#' reference system so distances can be evaluated in metres.
#'
#' \code{num_threads} uses OpenMP when available. On macOS, installing OpenMP
#' support with \code{brew install libomp} before installing the package may be
#' required for multi-threaded execution.
#'
#' @param data A matrix, data.frame, or \pkg{terra} \code{SpatRaster} containing
#' input RS data. Matrix and data.frame inputs must be organised as \code{x},
#' \code{y}, predicted RS variables, then observed RS variables. Raster inputs
#' must contain predicted RS layers followed by observed RS layers in the same
#' variable order;
#' @param samples A two-column matrix or data.frame of sample coordinates used
#' to extract predicted and observed RS values when \code{data} is a raster.
#' Ignored for matrix or data.frame inputs.
#' @param radius_km Numeric. Search radius, in kilometres, for deciding which
#' reference-sample pairs contribute to the density surface.
#' @param bin_width Numeric. Width of each predicted-distance and
#' observed-distance bin. Smaller values increase resolution but can produce a
#' sparse surface; larger values smooth more aggressively. The value is stored
#' as a \code{bin.width} attribute on the output.
#' @param bin_num Integer. Number of bins along each axis of the density surface.
#' The default is usually adequate; tuning \code{bin_width} is generally more
#' useful than changing \code{bin_num}.
#' @param drop_features Optional integer vector of RS variable positions to
#' exclude from the density calculation. Positions are 1-based within the RS
#' feature set, not within the full input column order. Use the same exclusion
#' in \code{\link{benchmark}} unless there is a deliberate reason not to.
#' @param num_threads Integer. Number of CPU threads to use. Values below 1 use
#' all available OpenMP threads.
#' @param filename Optional character. File path for writing the raw density
#' surface as a tab-delimited \file{.txt} file.
#'
#' @seealso \code{\link{normalise}}, and \code{\link{benchmark}}
#'
#' @return A \code{reference_density} object, which is also a matrix/array.
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' # Matrix inputs are x, y, predicted RS variables, then observed RS variables.
#' reference_data <- cbind(
#'     x = c(0, 0.4, 0.8, 1.2, 1.6, 2.0),
#'     y = c(0, 0.2, 0.8, 1.0, 1.4, 1.8),
#'     rs1 = c(0.10, 0.12, 0.25, 0.30, 0.42, 0.50),
#'     rs1 = c(0.11, 0.14, 0.22, 0.33, 0.40, 0.52)
#' )
#'
#' rd <- ref_density(
#'     reference_data,
#'     radius_km = 250,
#'     bin_width = 0.1,
#'     bin_num = 20,
#'     num_threads = 1
#' )
#'
#' rd_norm <- normalise(rd, trim_size = 10)
#' }
ref_density <- function(
        data,
        samples = NULL,
        radius_km = 1000,
        bin_width = 0.05,
        bin_num = 650,
        drop_features = NULL,
        num_threads = -1,
        filename = "") {

    if (radius_km <= 0) stop("radius_km must a postive non-zero number.")

    # check for data variables
    if (.is_mat(data)) {
        data_vals <- .check_mat(data)
        keep_features <- .keep_rs_features(drop_features, .num_rs_vars_mat(data_vals, "data"))
        data_vals <- .subset_rs_mat(data_vals, keep_features)
    } else if (.is_rast(data)) {
        # check terra is available
        .check_pkgs("terra")

        # check samples
        if (is.null(samples)) stop("For input 'data' as a raster file, 'sample' xy must be provided!")

        samples <- if (.is_mat(samples)) .check_mat(samples) else stop("'samples' must be a matrix or convertible to one.")

        if (ncol(samples) != 2) stop("'samples' must be a data.frame or matrix with exactly two columns of XY coordinates.")

        data <- .check_rast(data)
        if (terra::nlyr(data) %% 2) stop("Odd number of layers! The number of observed and prediced RS must be equal!")

        keep_features <- .keep_rs_features(drop_features, terra::nlyr(data) / 2L)
        data <- .subset_rs_rast(data, keep_features)

        # extract values
        data_ext <- terra::extract(data, samples, ID = FALSE)
        data_vals <- as.matrix(cbind(samples, data_ext))
    }  else {
        stop("The 'predicted' must be raster or a matrix, or convertiable object to these classes.")
    }

    # number of observed RS vars
    num_layers <- .num_rs_vars_mat(data_vals, "data")
    # id of obs and mod for saving
    mod_layers <- seq_len(num_layers) + 2
    obs_layers <- mod_layers + num_layers

    # get the correct columns for the C++ code
    samples_xy <- data_vals[, 1:2, drop = FALSE]
    modelled <- data_vals[, mod_layers, drop = FALSE]
    observed <- data_vals[, obs_layers, drop = FALSE]

    # some error checking
    if(any(dim(modelled) != dim(observed)))
        stop("Dimensions of RS and ENV datasets doesn't match!")

    if(nrow(modelled) != nrow(samples_xy))
        stop("Number of rows of rasters values and reference samples_xy doesn't match!")

    if(any(anyNA(modelled), anyNA(observed), anyNA(samples_xy)))
        stop("There's NA in the extracted observed, modelled or xy samples_xy!")

    if(any(colnames(modelled) != colnames(observed))) {
        warning("The names\\order of observed and modelled datasets doesn't match!\n")
        cat("Observed:", colnames(observed), "\n")
        cat("Modelled:", colnames(modelled), "\n")
    }

    # run reference density calculation in C++
    tryCatch(
        {
            out_table <- ref_density_cpp(
                rs_vals = observed,
                pr_vals = modelled,
                xy_vals = samples_xy,
                radius_km = radius_km,
                bin_width = bin_width,
                bin_num = bin_num,
                geographic = .is_lonlat(samples_xy),
                num_threads = num_threads
            )
        },
        error = function(cond) {
            stop("Reference density calculation failed!\n", cond)
        }
    )

    # write down the reference density surface
    if (nchar(filename) > 0) {
        filename <- ifelse(grepl(".txt", filename), filename, paste0(filename, ".txt"))
        tryCatch(
            {
                utils::write.table(
                    out_table,
                    file = filename,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE
                )
            },
            error = function(cond) {
                stop("Writing reference density file failed!\n", cond)
            }
        )
    }

    class(out_table) <- c("reference_density", "matrix", "array")
    attr(out_table, "bin.width") <- bin_width

    return(out_table)
}


#' @export
#' @method print reference_density
print.reference_density <- function(x, ...) {
    print(class(x), ...)
    i <- which(names(attributes(x)) == "class")
    print(attributes(x)[-i])
}

#' @export
#' @method plot reference_density
plot.reference_density <- function(x, ...) {
    is_norm <- "offset" %in% names(attributes(x))
    if (is_norm) {
        x <- apply(t(x), 2, rev)
        message("For aesthetic, the normalised reference density plot is reversed and transposed.")
    }
    terra::plot(
        terra::rast(.check_mat(unclass(x), name = "x")), col = palettes(150, "ref_density"), ...
    )
}
