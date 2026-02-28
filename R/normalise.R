#' Clean and normalise HCAS reference density
#'
#' This process includes trimming the reference density surface to remove noise and
#' normalising its values.
#'
#' @param x An HCAS \code{reference_density} object or a matrix representing the
#' reference density surface (see \code{\link{ref_density}}).
#' @param bin_width Numeric. Specifies the bin width of the reference density. If
#' \code{x} is a \strong{reference_density} object, this value can be read from
#' its attributes and may be left \code{NULL}. The bin width must be consistent
#' between the reference density creation and the benchmarking step to ensure
#' condition is accurately calculated.
#' @param trim_size Integer. Defines the number of rows and columns in the trimmed reference density. The default
#' is 400, and it is generally advisable to retain this default setting.
#' @param offset Integer. Specifies the number of reference density bins to ignore during normalization. This value
#' will be stored as an attribute in the output object.
#' @param legacy Logical. Whether to use the legacy C++ code for normalisation (for backward
#' compatibility) or the modern R version (default). The modern version solves the edge effect
#' issue without any speed compromise.
#' @param filename Char (optional). The output file name for the .text file.
#'
#' @seealso \code{\link{ref_density}}, and \code{\link{benchmark}}
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
normalise <- function(
        x,
        bin_width = NULL,
        trim_size = 400,
        offset = 0,
        legacy = FALSE,
        filename = "") {

    # check bin_width and get it from reference density object
    if (methods::is(x, "reference_density")) {
        if (is.null(bin_width)) {
            bin_width <- attributes(x)$bin.width
        } else {
            if (bin_width != attributes(x)$bin.width) {
                warning("The supplied 'bin_width` is different from the arrtibute(x)$bin.width from the input.")
            }
        }
    }

    # force offset to be above zero
    offset <- max(0, offset)

    # remove the 0,0 point; self overlaps
    nr <- nrow(x)
    x[nr, 1] <- 0

    # normalise the reference density surface, C++ or R
    out <- if (legacy) {
        norm_cpp(x, trim_size = trim_size, offset = offset)
    } else {
        norm_r(x, trim_size = trim_size, offset = offset)
    }

    # write down the reference density surface
    if (nchar(filename) > 0) {
        filename <- ifelse(grepl(".txt", filename), filename, paste0(filename, ".txt"))
        tryCatch(
            {
                utils::write.table(
                    out,
                    file = filename,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE
                )
            },
            error = function(cond) {
                stop("Writing reference density file failed!\n", cond)
            }
        )
    }

    class(out) <- c("reference_density", "matrix", "array")
    if (is.null(bin_width)) {
        attr(out, "bin.width") <- NA
    } else {
        attr(out, "bin.width") <- bin_width
    }
    attr(out, "offset") <- offset

    return(
        out
    )
}

# the new normalisation function that doesn't remove the edges
norm_r <- function(x, trim_size = 400, offset = 0) {
    nr <- nrow(x)
    nc <- ncol(x)

    r <- terra::rast(x)
    # create a Gaussian filter and normalise it to add up to 1
    d <- stats::dnorm(-2:2, 0, 1)
    w <- outer(d, d)
    w <- w / sum(w)
    rr <- terra::focal(r, w = w, fun = sum, na.rm = TRUE)
    mat <- terra::as.matrix(rr, wide = TRUE)

    # trim the reference density and apply offset
    upper <- nr - offset
    lower <- min(nr - trim_size + 1 - offset, upper)
    left <- 1 + offset
    right <- min(trim_size + offset, nc)
    mat <- mat[lower:upper, left:right]

    # normalise the columns and then reverse them
    mat <- apply(mat, 2, FUN = function(x) rev(x / sum(x)))

    return(
        t(mat)
    )
}

