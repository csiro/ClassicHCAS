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

    # normalise the histo C++
    out <- norm_cpp(
        x = x,
        trim_size = trim_size,
        offset = offset
    )

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

