#' Reference density surface
#'
#' The HCAS reference density surface calculation is based on pair-point densities.
#' It is an integral part of HCAS, designed to
#' learn the expected observed remote sensing (RS) values from the predicted RS
#' values across a wide range of reference sites. Note that the reference sites
#' used for \code{ref_density} do not need to be the same as those used in the
#' \code{\link{benchmark}} function. See more in details.
#'
#' Ensure that the order of remote sensing variables is consistent between predicted and observed inputs
#' (for both raster and matrix formats). The RS variable values must be centered and scaled
#' prior to prediction. Failure to do so may result in variables with larger ranges having
#' disproportionate influence in the multi-dimensional distance calculations.
#'
#' This function uses an integer-based distance checks for fast radius searches on either
#' geographic or projected coordinates. In geographic mode, coordinates are
#' stored in micro-degrees (degree * 1000_000) and the distance is approximated by:
#'
#'     distance² ≈ (dlat)² + (dlon × cos(lat₁))²
#'
#' where cos(lat₁) is derived from the query latitude. This avoids floating-
#' point overhead and provides substantial performance gains but introduces
#' distortion at larger distances. For applications requiring higher accuracy,
#' especially beyond regional scales (more than several 100s of kilometers in \code{radius_km}),
#' use a projected coordinate system so distances in meters can be evaluated directly.
#'
#' Ensure that \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system
#' to take advantage of parallel processing and accelerate computations. While most systems
#' include OpenMP by default, you may need to load the appropriate module if you're using an HPC
#' system.
#'
#' \strong{Note for macOS users:} Install OpenMP via Homebrew with \code{brew install libomp}
#' before installing this package.
#'
#' @param data A matrix, SpatRaster (from the \pkg{terra} package), or data.frame containing the input data.
#' The data \strong{must} be organised in the following order: \strong{x}, \strong{y}, \strong{predicted-RS},
#' \strong{observed-RS} variables. If using a SpatRaster, the \strong{x} and \strong{y} are not required.
#' If using a matrix or data.frame, ensure that the variables are in the
#' correct order. For more see the details section.
#' @param samples A matrix or data.frame containing the x and y coordinates (longitude and latitude)
#' of the reference points used for the observed and predicted RS value extraction if \code{data} is a
#' raster object. If \code{data} argument is matrix, this will be ignored.
#' @param radius_km Numeric. Specifies the search radius in kilometers for considering reference samples
#' when creating the reference density surface. See details section for more information on distance calculation.
#' @param bin_width Numeric. Specifies the bin width of the reference density surface. Finding the optimal bin width
#' may require some experimentation to achieve the best results. The bin width is added as an attribute
#' to the output object, ensuring consistency and accuracy in subsequent benchmarking steps.
#' @param bin_num Integer. Specifies the number of bins for the reference density surface. It is generally recommended
#' to use the default value of 650. Adjusting \code{bin_width} is often more effective than changing
#' \code{bin_num}.
#' @param drop_features Integer vector. Completely remove the RS variable from the reference density calculation. For
#' consistency, it is recommended to exclude the same variables later in the benchmarking step; unless
#' you have a specific reason not to.
#' @param num_threads Integer. Specifies the number of CPU threads to be used for processing. A value
#' below 1 indicates that all available threads will be utilized (default). Refer to the details section for
#' more information.
#' @param filename Character (optional). The name of the output file for saving the results as a .txt file.
#'
#' @seealso \code{\link{normalise}}, and \code{\link{benchmark}}
#'
#' @return A \code{reference_density} object (also matrix, array)
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#'
#'
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
    } else if (.is_rast(data)) {
        # check terra is available
        .check_pkgs("terra")

        # check samples
        if (is.null(samples)) stop("For input 'data' as a raster file, 'sample' xy must be provided!")

        samples <- if (.is_mat(samples)) .check_mat(samples) else stop("'samples' must be a matrix or convertible to one.")

        if (ncol(samples) != 2) stop("'samples' must be a data.frame or matrix with exactly two columns of XY coordinates.")

        data <- .check_rast(data)
        if (terra::nlyr(data) %% 2) stop("Odd number of layers! The number of observed and prediced RS must be equal!")

        # extract values
        data_ext <- terra::extract(data, samples, ID = FALSE)
        data_vals <- as.matrix(cbind(samples, data_ext))
    }  else {
        stop("The 'predicted' must be raster or a matrix, or convertiable object to these classes.")
    }


    # number of observed RS vars
    num_layers <- (ncol(data_vals) - 2) / 2
    # id of obs and mod for saving
    mod_layers <- seq_len(num_layers) + 2
    obs_layers <- mod_layers + num_layers

    # get the correct columns for the C++ code
    samples_xy <- data_vals[, 1:2]
    modelled <- data_vals[, mod_layers]
    observed <- data_vals[, obs_layers]

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

    # drop features from calculation if requested
    if (length(drop_features)) {
        observed[, drop_features] <- 0
        modelled[, drop_features] <- 0
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
        terra::rast(x), col = ref_density_color(150), ...
    )
}

