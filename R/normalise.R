#' Clean and normalise HCAS histogram
#'
#' This process includes trimming the histogram to remove noise and normalising its values.
#'
#' @param x An HCAS histo object or a matrix representing the HCAS histogram (see \code{\link{histogram}}).
#' @param bin_width Numeric. Specifies the bin width of the histogram. If \code{x} is a
#' \strong{histo} object, this value can be read from its attributes and may be left \code{NULL}. The bin width
#' is primarily used to ensure consistency in the benchmarking step by accurately tracking and applying the
#' bin width used during histogram creation.
#' @param trim_size Integer. Defines the number of rows and columns in the trimmed histogram. The default
#' is 400, and it is generally advisable to retain this default setting.
#' @param offset Integer. Specifies the number of histogram bins to ignore during normalization. This value
#' will be stored as an attribute in the output object.
#' @param legacy Logical. Whether to use the legacy C++ code for normalisation (for backward
#' compatibility) or the modern R version (default). The modern version solves the edge effect issue.
#' @param filename Char (optional). The output file name for the .text file.
#'
#' @seealso \code{\link{histogram}}, and \code{\link{benchmark}}
#'
#' @return A histo object (matrix, array)
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

    # check bin_width and get it from histo object
    if (methods::is(x, "histo")) {
        if (is.null(bin_width)) {
            bin_width <- attributes(x)$bin.width
        } else {
            if (bin_width != attributes(x)$bin.width) {
                warning("The supplied 'bin_width` is different from the arrtibute(x)$bin.width from the input.")
            }
        }
    }

    # remove the 0,0 point
    nr <- nrow(x)
    x[nr, 1] <- 0

    # normalise the histo C++
    if (legacy) {
        out <- norm_cpp(
            x = x,
            trim_size = trim_size,
            offset = offset
        )
    } else {
        out <- norm_r(
            x = x,
            trim_size = trim_size,
            offset = offset
        )
    }

    # write down the histogram
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
                stop("Writing histogram file failed!\n", cond)
            }
        )
    }

    class(out) <- c("matrix", "array", "histo")
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
    x[nr, 1] <- 0
    r <- terra::rast(x)
    # create a Gaussian filter and normalise it to add up to 1
    d <- stats::dnorm(-2:2, 0, 1)
    w <- outer(d, d)
    w <- w / sum(w)
    rr <- terra::focal(r, w = w, fun = sum, na.rm = TRUE)
    mat <- as.matrix(rr, wide = TRUE)

    # trim the histogram and apply offset
    lower <- nr - trim_size + 1 - offset
    upper <- nr - offset
    left <- 1 + offset
    right <- trim_size + offset
    mat <- mat[lower:upper, left:right]

    # normalise the columns and then reverse them
    mat <- apply(mat, 2, FUN = function(x) rev(x / sum(x)))

    return(
        t(mat)
    )
}

