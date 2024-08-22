#' Clean and normalise HCAS histogram
#'
#' This process includes trimming the histogram to remove noise and normalising the values.
#'
#' @param x An HCAS histo object or a matrix representing the HCAS histogram (see \code{\link{histogram}}).
#' @param bin_width Integer (optional). Specifies the bin size of the histogram. If \code{x} is a
#' \strong{histo} object, this value can be read from its attributes and may be left \code{NULL}. The bin size
#' is primarily used to ensure consistency in the benchmarking step by accurately tracking and applying the
#' bin width used during histogram creation.
#' @param trim_size Integer. Defines the number of rows and columns in the trimmed histogram. The default
#' is 400, and it is generally advisable to retain this default setting.
#' @param offset Integer. Specifies the number of histogram bins to ignore during normalization. This value
#' will be stored as an attribute in the output object.
#'
#' @return A histo object (matrix, array)
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
clean <- function(x, bin_width = NULL, trim_size = 400, offset = 0) {
    # check bin_width and get it from histo object
    if (is(x, "histo")) {
        if (is.null(bin_width)) {
            bin_width <- attributes(x)$bin.width
        } else {
            if (bin_width != attributes(x)$bin.width) {
                warning("The supplied 'bin_width` is different from the arrtibute(x)$bin.width from the input.")
            }
        }
    }

    out <- clean_cpp(
        x = x,
        trim_size = trim_size,
        offset = offset
    )

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

