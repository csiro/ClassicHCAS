#' Clean and normalise an HCAS reference density surface
#'
#' Trims and normalises the raw reference density surface returned by
#' \code{\link{ref_density}} so it can be used as a probability surface in
#' \code{\link{benchmark}}.
#'
#' @details
#' The raw reference density surface is a two-dimensional surface of predicted
#' and observed RS distances among reference samples. Before benchmarking, the
#' surface is smoothed, trimmed to remove noisy outer bins, and normalised with
#' respect to predicted-distance bins. This makes each predicted-distance slice
#' comparable when \code{\link{benchmark}} asks how probable an observed
#' departure is for a target location.
#'
#' \code{trim_size} controls the dimensions of the retained square surface. It
#' should be smaller than the number of bins used in
#' \code{\link{ref_density}}. The default is chosen for the standard HCAS
#' workflow; smaller examples or exploratory analyses can use smaller values.
#'
#' \code{offset} removes bins nearest the origin before normalisation. It is
#' useful when the near-zero distance cells contain self-overlap or other
#' artefacts. The value is stored on the returned \code{reference_density}
#' object and must be kept consistent in \code{\link{benchmark}}.
#'
#' @param x An HCAS \code{reference_density} object or a matrix representing the
#' raw reference density surface created by \code{\link{ref_density}}.
#' @param bin_width Numeric. Bin width used to create the reference density. If
#' \code{x} is a \code{reference_density} object, the value is read from its
#' \code{bin.width} attribute when \code{bin_width = NULL}. The value must match
#' the density surface used during benchmarking.
#' @param trim_size Integer. Number of rows and columns to keep in the trimmed
#' reference density surface.
#' @param offset Integer. Number of near-origin bins to ignore during
#' normalisation. Stored as an attribute on the output.
#' @param legacy Logical. If \code{TRUE}, use the legacy C++ normalisation code
#' for backward compatibility. The default R implementation avoids the previous
#' edge effect while retaining similar speed for typical use.
#' @param filename Optional character. File path for writing the normalised
#' surface as a tab-delimited \file{.txt} file.
#'
#' @seealso \code{\link{ref_density}}, and \code{\link{benchmark}}
#'
#' @return A \code{reference_density} object, which is also a matrix/array.
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' raw <- matrix(rexp(30 * 30), nrow = 30)
#' class(raw) <- c("reference_density", "matrix", "array")
#' attr(raw, "bin.width") <- 0.1
#'
#' norm <- normalise(raw, trim_size = 15)
#' attr(norm, "bin.width")
#' attr(norm, "offset")
#' }
normalise <- function(
        x,
        bin_width = NULL,
        trim_size = 400,
        offset = 0,
        legacy = FALSE,
        filename = "") {

    if (!.is_mat(x)) {
        stop("'x' must be a matrix or a 'reference_density' object.")
    }

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

    # drop the custom class so terra and Rcpp see a plain matrix
    x <- .check_mat(unclass(x), name = "x")

    # force offset to be above zero
    offset <- max(0, offset)

    # remove the 0,0 point; self overlaps
    nr <- nrow(x)
    x[nr, 1] <- 0

    # normalise the reference density surface, C++ or R
    out <- if (legacy) {
        norm_cpp(x, trim_size = trim_size, offset = offset)
    } else {
        .norm_r(x, trim_size = trim_size, offset = offset)
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
.norm_r <- function(x, trim_size = 400, offset = 0) {
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
