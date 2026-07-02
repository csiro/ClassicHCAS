#' Benchmark target locations against HCAS reference condition
#'
#' Estimates raw HCAS habitat condition by comparing target locations with
#' nearby, environmentally similar reference samples and a normalised reference
#' density surface.
#'
#' @details
#' In HCAS, predicted RS variables represent the expected signal under reference
#' condition and observed RS variables represent the actual Earth observation
#' signal. \code{benchmark()} asks whether the observed departure from expected
#' reference condition is typical of high-integrity reference ecosystems. The
#' output is an unscaled relative condition value; use \code{\link{calibrate}} to
#' map it to a 0-1 condition scale.
#'
#' \strong{Experimental temporal mode:} the mechanism controlled by
#' \code{temporal_sigma} is exploratory, is off by default, and should not be
#' treated as equivalent to the standard non-temporal benchmark without
#' independent validation. Standard benchmarking remains non-temporal unless a
#' finite \code{temporal_sigma} is supplied.
#'
#' The function uses a two-stage reference-sample selection process for each
#' target location:
#' \enumerate{
#'   \item Candidate benchmark samples are restricted to those within
#'   \code{radius_km}.
#'   \item From those candidates, up to \code{k1} samples with the smallest
#'   predicted RS distance are retained. If \code{xy_penalty > 0}, scaled
#'   geographic coordinates are included in this distance so distant samples are
#'   penalised even when they are spectrally similar.
#'   \item The target and retained samples are queried against
#'   \code{ref_density} using predicted-distance and observed-distance bins.
#'   If the experimental temporal mode is enabled, each retained site is queried
#'   once per reference year, and a Gaussian year weight selects and weights the
#'   site's most relevant year.
#'   \item Up to \code{k2} samples with the highest reference-density
#'   probability are retained for condition estimation.
#' }
#'
#' The retained probability values are combined using the distance kernel
#' selected by \code{kernel}. The default Gaussian kernel is
#' \code{exp(-(distance / lambda)^2)}. The optional Cauchy kernel is
#' \code{1 / (1 + (distance / lambda)^2)}, the standard Cauchy shape
#' normalised to weight one at zero.
#'
#' By default, \code{boost = 10} replaces the LDC blend with a boosted weighted
#' mean. The kernel weight of the retained site with the highest unweighted
#' probability is multiplied by \code{boost}, and condition is the weighted mean
#' using that adjusted weight. In this mode, \code{confidence} does not affect
#' condition. Set \code{boost = NULL} or \code{boost = NA} to use the original
#' LDC blend, where \code{confidence} controls how strongly the raw condition
#' value relies on the maximum unweighted probability contribution compared
#' with the distance-weighted mean probability. A value of \code{boost = 1}
#' gives the ordinary kernel-weighted mean. If
#' \code{make_su = TRUE}, the result also includes \code{su}, the log of the
#' original, unboosted total distance-weight sum, which is a support diagnostic
#' rather than a calibrated confidence interval.
#'
#' Matrix and data.frame inputs must be ordered as \code{x}, \code{y}, predicted
#' RS variables, then observed RS variables. Raster inputs must contain predicted
#' RS layers followed by observed RS layers in the same variable order. The
#' sample matrix can either contain only \code{x} and \code{y} coordinates, in
#' which case values are extracted from a raster \code{data} object, or the full
#' \code{x}, \code{y}, predicted RS, observed RS table. Predicted and observed
#' variables should be centred and scaled consistently before benchmarking.
#'
#' The most influential tuning parameters are usually \code{radius_km},
#' \code{xy_penalty}, \code{k1}, \code{k2}, \code{lambda}, and
#' \code{boost}. Defaults were chosen for Australian HCAS applications and
#' should be assessed before use in other regions, data products, or ecological
#' contexts. When field condition data are unavailable, tuning can be guided by
#' whether scores discriminate among independent land-use or disturbance classes
#' in the expected order.
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
#' @inheritParams ref_density
#' @param samples Benchmark/reference sample data. Normally this is a matrix or
#' data.frame containing \code{x}, \code{y}, predicted RS variables, then
#' observed RS variables in the same order as \code{data}. For raster
#' \code{data}, it can instead contain only two coordinate columns, in which
#' case raster values are extracted. When \code{temporal_sigma} is specified,
#' supply a named list of full sample matrices, one per year, with numeric year names.
#' All matrices must contain the same sites, row order, XY coordinates,
#' predicted values, and feature order; only observed values may vary.
#' @param ref_density A normalised \code{reference_density} object or matrix
#' produced by \code{\link{normalise}}.
#' @param xy_stats Numeric vector of length four used to scale coordinates when
#' \code{xy_penalty > 0}: \code{mean(x)}, \code{mean(y)}, \code{sd(x)},
#' \code{sd(y)}. Use the same values across tiles to keep tiled benchmarking
#' consistent.
#' @param xy_penalty Numeric. Weight applied to scaled coordinates when selecting
#' the \code{k1} most similar benchmark samples. \code{0} disables the
#' spatial penalty.
#' @param radius_km Numeric. Search radius, in kilometres, for candidate
#' benchmark samples.
#' @param k1 Integer. First-stage filter size: the number of nearest samples to
#' retain after the predicted RS distance search.
#' @param k2 Integer. Second-stage filter size: the number of high-probability
#' samples to retain from the reference density query. Must be less than or
#' equal to \code{k1}.
#' @param bin_width Numeric. Bin width used to create and normalise
#' \code{ref_density}. If \code{ref_density} is a \code{reference_density}
#' object, this value is read from its \code{bin.width} attribute when
#' \code{bin_width = NULL}.
#' @param interpolate Logical. If \code{TRUE}, bilinearly interpolates the
#' reference density surface before benchmarking for smoother lookup.
#' @param offset Integer. Number of reference-density bins ignored during
#' normalisation. If \code{ref_density} is a \code{reference_density} object,
#' this value is read from its \code{offset} attribute when \code{offset = NULL}.
#' @param confidence Numeric between 0 and 1. Weight given to the selected
#' maximum probability component relative to the distance-weighted mean
#' probability when computing raw condition. Ignored when \code{boost} is not
#' \code{NULL} or \code{NA}; the default is \code{boost = 10}.
#' @param boost \code{NULL}, \code{NA}, or one positive finite numeric factor.
#' The default \code{10} multiplies the kernel weight of the
#' highest-probability retained site by ten and returns the resulting weighted
#' mean instead of the LDC blend. \code{confidence} is ignored in this mode.
#' Use \code{NULL} or \code{NA} for the unboosted LDC blend.
#' @param lambda Positive numeric. Distance-scale bandwidth for the selected
#' \code{kernel}. Both kernels treat \code{lambda} in predicted RS L1 distance
#' units: the Gaussian kernel uses \code{exp(-(distance / lambda)^2)} and the
#' Cauchy kernel uses \code{1 / (1 + (distance / lambda)^2)}.
#' @param exclude_slef Logical. If \code{TRUE}, exclude samples whose predicted
#' RS distance is less than one bin width, preventing a benchmark point from
#' assessing itself. The argument name preserves the existing API spelling.
#' @param drop_features Optional integer vector of RS variable positions to
#' exclude from benchmarking. Positions are 1-based within the RS feature set,
#' not within the full input column order. Use the same exclusion used in
#' \code{\link{ref_density}} unless there is a deliberate reason not to.
#' @param assessment_year Numeric. Year being assessed and the centre of the
#' experimental Gaussian temporal kernel. Required when
#' \code{temporal_sigma} is specified.
#' @param temporal_sigma Positive numeric or \code{NULL}. Standard deviation of
#' the experimental Gaussian temporal kernel, in the same units as the
#' temporal sample names. Supplying a value enables an exploratory temporal
#' mode that has not been validated as a drop-in replacement for standard
#' non-temporal benchmarking; \code{NULL} or \code{NA} disables it. When
#' enabled, the selected reference-density probability is multiplied by the
#' selected year's temporal weight before \code{k2} selection and condition
#' estimation.
#' @param make_su Logical. If \code{TRUE}, return both raw condition and
#' \code{su}, the log of the total distance-weight sum.
#' @param kernel Character. Distance kernel applied to retained-reference
#' predicted RS L1 distances: \code{"Gaussian"} (default) or \code{"Cauchy"}.
#' Lower-case \code{"gaussian"} and \code{"cauchy"} are also accepted.
#' @param ... Additional arguments passed to \code{\link[terra]{interpolate}}
#' when benchmarking raster outputs, such as \code{filename}, \code{overwrite},
#' or \code{wopt}.
#'
#' @seealso \code{\link{ref_density}}, \code{\link{normalise}}, and \code{\link{calibrate}}
#'
#' @return A matrix or \pkg{terra} \code{SpatRaster}, depending on the inputs.
#' @export
#'
#' @examples
#' \donttest{
#' library(ClassicHCAS)
#'
#' target_data <- cbind(
#'     x = c(0.1, 0.9),
#'     y = c(0.1, 0.9),
#'     rs1 = c(0.12, 0.42),
#'     rs1 = c(0.13, 0.50)
#' )
#'
#' sample_data <- cbind(
#'     x = c(0.0, 0.4, 0.8, 1.2),
#'     y = c(0.0, 0.3, 0.8, 1.1),
#'     rs1 = c(0.10, 0.20, 0.40, 0.55),
#'     rs1 = c(0.11, 0.18, 0.43, 0.58)
#' )
#'
#' ref <- matrix(1, nrow = 20, ncol = 20)
#' class(ref) <- c("reference_density", "matrix", "array")
#' attr(ref, "bin.width") <- 0.1
#' attr(ref, "offset") <- 0
#'
#' benchmark(
#'     target_data,
#'     samples = sample_data,
#'     ref_density = ref,
#'     radius_km = 200,
#'     k1 = 3,
#'     k2 = 2,
#'     interpolate = FALSE,
#'     exclude_slef = FALSE,
#'     num_threads = 1
#' )
#' }
benchmark <- function(
        data,
        samples,
        ref_density,
        xy_stats = c(0, 0, 1, 1),
        xy_penalty = 0.0,
        radius_km = 200,
        k1 = 50,
        k2 = 20,
        bin_width = NULL,
        interpolate = TRUE,
        offset = 0,
        kernel = c("Gaussian", "Cauchy"),
        lambda = 1.0,
        confidence = 0.5,
        exclude_slef = TRUE,
        drop_features = NULL,
        assessment_year = NULL,
        temporal_sigma = NULL,
        make_su = FALSE,
        num_threads = -1,
        boost = 10,
        ...) {

    kernel <- .check_kernel(kernel)
    dots <- list(...)
    legacy_k <- intersect(names(dots), c("k_pred", "k_obs"))
    if (length(legacy_k)) {
        stop("'k_pred' and 'k_obs' were renamed to 'k1' and 'k2'.")
    }
    if ("temporal_correct" %in% names(dots)) {
        stop(
            "'temporal_correct' has been removed; supply 'temporal_sigma' ",
            "to enable the experimental temporal mode."
        )
    }
    if ("temporal_weighted" %in% names(dots)) {
        stop(
            "'temporal_weighted' has been removed; temporal weights are ",
            "always applied when 'temporal_sigma' enables the experimental ",
            "temporal mode."
        )
    }
    if ("weighted_max" %in% names(dots)) {
        stop("'weighted_max' is not an argument to benchmark().")
    }

    # The second-stage filter cannot retain more samples than the first stage.
    if (k1 < k2) stop("'k2' must be less than or equal to 'k1'.")
    boost <- .check_boost(boost)

    temporal_sigma_missing <- is.null(temporal_sigma) ||
        (length(temporal_sigma) == 1L && isTRUE(is.na(temporal_sigma)))
    use_temporal_correction <- !temporal_sigma_missing

    temporal <- NULL
    temporal_weights <- NULL
    if (use_temporal_correction) {
        temporal <- .prepare_temporal_samples(
            samples = samples,
            assessment_year = assessment_year,
            temporal_sigma = temporal_sigma,
            drop_features = drop_features
        )
        samples <- temporal$samples
        temporal_weights <- temporal$weights
    } else {
        samples <- if (.is_mat(samples)) .check_mat(samples) else stop("'samples' must be a matrix or convertible to one.")
    }

    # check reference density
    ref_density <- if (.is_mat(ref_density)) .check_mat(ref_density) else stop("'ref_density' must be a matrix or convertible to one.")
    if (nrow(ref_density) != ncol(ref_density)) warning("Reference density dimensions are not equal!\n")

    if (methods::is(ref_density, "reference_density")) {
        # check for reference density bin_width consistency
        if (is.null(bin_width)) {
            bin_width <- attributes(ref_density)$bin.width
        } else {
            if (bin_width != attributes(ref_density)$bin.width) {
                warning("Provided 'bin_width' differs from reference density attribute.")
            }
        }
        # check for reference density offset consistency
        if (is.null(offset)) {
            offset <- attributes(ref_density)$offset
        } else {
            if (offset != attributes(ref_density)$offset) {
                warning("Provided 'offset' differs from reference density attribute.")
            }
        }
    }

    # interpolate reference density
    if (interpolate) {
        ref_density <- terra::as.matrix(
            terra::disagg(
                terra::rast(.check_mat(unclass(ref_density), name = "ref_density")),
                fact = 2,
                method = "bilinear"
            ),
            wide = TRUE
        )
        # update the binwidth and offset
        bin_width <- bin_width / 2
        offset <- offset * 2
    }
    # get the bin number after interpolation
    bin_num <- min(dim(ref_density))

    if (.is_mat(data)) {
        # check and convert to matrix
        data <- .check_mat(data)

        if (use_temporal_correction) {
            if (.num_rs_vars_mat(data, "data") != temporal$n_vars) {
                stop("Temporal samples and 'data' must contain the same RS feature count.")
            }
            data <- .subset_rs_mat(data, temporal$keep_features)
        } else {
            if (ncol(samples) != ncol(data)) {
                stop("Samples must include all raster values (matching column count with 'data').")
            }

            keep_features <- .keep_rs_features(drop_features, .num_rs_vars_mat(samples, "samples"))
            samples <- .subset_rs_mat(samples, keep_features)
            data <- .subset_rs_mat(data, keep_features)
        }

        tryCatch(
            {
                output <- .benchmarking(
                    model = list(),
                    newdata = data, # rast_stack arg
                    sample_vals = samples,
                    ref_density = ref_density,
                    xy_stats = xy_stats,
                    xy_penalty = xy_penalty,
                    radius_km = radius_km,
                    geographic = .is_lonlat(data),
                    bin_width = bin_width,
                    bin_num = bin_num,
                    offset = offset,
                    k_env = k1,
                    k_rs = k2,
                    confidence = confidence,
                    boost = boost,
                    lambda = lambda,
                    exclude_slef = exclude_slef,
                    temporal_weights = temporal_weights,
                    make_su = make_su,
                    num_threads = num_threads,
                    kernel = kernel
                )
            },
            error = function(cond) {
                stop("HCAS benchmarking C++ function failed!\n", cond)
            }
        )
    } else if (.is_rast(data)) {
        # check terra is available
        .check_pkgs("terra")
        # check and convert to SpatRaster object
        data <- .check_rast(data)

        if (terra::nlyr(data) %% 2L) {
            stop("'data' must contain matching predicted and observed RS layers.")
        }

        if (use_temporal_correction) {
            if (terra::nlyr(data) / 2L != temporal$n_vars) {
                stop("Temporal samples and 'data' must contain the same RS feature count.")
            }
            data <- .subset_rs_rast(data, temporal$keep_features)
        } else {
            # sample extraction if needed
            if (ncol(samples) == 2) {
                cat("Extracting sample values...\n")
                samples <- cbind(samples, as.matrix(terra::extract(data, samples, ID = FALSE)))
            } else if ((ncol(samples) - 2) != terra::nlyr(data)) {
                stop("Sample feature count does not match number of raster layers.")
            }

            keep_features <- .keep_rs_features(drop_features, terra::nlyr(data) / 2L)
            samples <- .subset_rs_mat(samples, keep_features)
            data <- .subset_rs_rast(data, keep_features)
        }

        tryCatch(
            {
                output <- terra::interpolate(
                    object = data,
                    model = list(),
                    fun = .benchmarking,
                    sample_vals = samples,
                    ref_density = ref_density,
                    xy_stats = xy_stats,
                    xy_penalty = xy_penalty,
                    radius_km = radius_km,
                    geographic = .is_lonlat(data),
                    bin_width = bin_width,
                    bin_num = bin_num,
                    offset = offset,
                    k_env = k1,
                    k_rs = k2,
                    confidence = confidence,
                    boost = boost,
                    lambda = lambda,
                    exclude_slef = exclude_slef,
                    temporal_weights = temporal_weights,
                    make_su = make_su,
                    num_threads = num_threads,
                    kernel = kernel,
                    ...
                )
            },
            error = function(cond) {
                stop("HCAS benchmarking C++ function failed!\n", cond)
            }
        )

    } else {
        stop("The 'data' must be raster or a matrix, or convertiable object to these classes.")
    }

    return(
        output
    )
}


# a function to handling predicting with terra
.benchmarking <- function(model, newdata, make_su, ...){
    nr <- nrow(newdata)
    nc <- make_su + 1
    col_names <- c("condition", "su")[1:nc]

    dat <- as.matrix(newdata)

    tryCatch(
        {
            hcas_cond <- bench_cpp(
                raster_vals = dat,
                make_su = make_su,
                ...
            )
        },
        error = function(cond) {
            message("Benchmarking C++ function failed. Returning -0.02 for all cells.")
            # return error values -0.02
            return(
                matrix(-0.02, nrow = nr, ncol = nc, dimnames = list(NULL, col_names))
            )
        }
    )

    colnames(hcas_cond) <- col_names

    return(hcas_cond)
}
