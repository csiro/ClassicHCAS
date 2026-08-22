#' Assess reference-site use during HCAS benchmarking
#'
#' Records how often and how strongly each reference site is used across the
#' three stages of the HCAS benchmarking selection process.
#'
#' @details
#' \code{reference_use()} runs the same non-temporal reference-site selection
#' stages used by \code{\link{benchmark}}, but returns one row per reference
#' site instead of habitat-condition values. It is implemented separately from
#' the benchmarking routine so normal benchmarking does not incur diagnostic
#' bookkeeping costs.
#'
#' The returned columns represent:
#' \describe{
#'   \item{\code{"predicted"}}{Count selection among the \code{k1} nearest
#'   reference sites in predicted feature space after the radius restriction
#'   and optional XY penalty.}
#'   \item{\code{"density"}}{Count retention among the \code{k2} sites with
#'   the highest reference-density probability.}
#'   \item{\code{"condition"}}{Attribute each retained site's use according to
#'   its contribution to the condition estimator. The
#'   distance-weighted mean component uses normalised weights from the selected
#'   \code{kernel}. With \code{boost} (the default is \code{k2}), attribution
#'   uses boosted normalised kernel weights. Set \code{boost = NULL} or
#'   \code{boost = NA} to use the LDC attribution, where the LDC component is
#'   assigned to the site with the selected maximum probability contribution,
#'   controlled by \code{weighted_max}; ties split that component equally.
#'   Attribution weights sum to one for each evaluated target location.}
#' }
#'
#' The experimental temporal mode is not currently supported. Supply a single
#' reference sample matrix containing \code{x}, \code{y}, predicted RS
#' variables, then observed RS variables.
#'
#' Raster inputs are processed in blocks and accumulated into the reference-site
#' table, so the complete raster is not loaded into memory.
#'
#' @inheritParams benchmark
#' @param data A matrix, data.frame, or \pkg{terra} \code{SpatRaster}
#' containing target RS data. Matrix and data.frame inputs must be organised as
#' \code{x}, \code{y}, predicted RS variables, then observed RS variables.
#' Raster inputs must contain predicted RS layers followed by observed RS
#' layers in the same variable order.
#' @param samples A matrix or data.frame containing reference sites as
#' \code{x}, \code{y}, predicted RS variables, then observed RS variables in
#' the same order as \code{data}. For raster \code{data}, this can instead be a
#' two-column coordinate table, in which case values are extracted from the
#' raster. Temporal sample lists are not supported.
#' @param weighted_max Logical. If \code{FALSE} (the default), the LDC
#' attribution component uses the maximum unweighted probability. If
#' \code{TRUE}, it uses the maximum probability-times-distance-weight
#' contribution. This option affects attribution only when \code{boost = NULL}
#' or \code{boost = NA}.
#'
#' @return A data.frame with one row per reference site and four columns:
#' \itemize{
#'   \item \code{id}: 1-based row index in \code{samples}.
#'   \item \code{predicted}: number of selections among \code{k1}.
#'   \item \code{density}: number of retentions among \code{k2}.
#'   \item \code{condition}: summed distance-kernel/LDC attribution weight.
#' }
#'
#' @seealso \code{\link{benchmark}}
#' @export
#'
#' @examples
#' target <- matrix(
#'     c(
#'         0, 0, 0.1, 0.1,
#'         1, 1, 0.8, 0.8
#'     ),
#'     ncol = 4,
#'     byrow = TRUE
#' )
#' samples <- matrix(
#'     c(
#'         0, 0, 0.1, 0.1,
#'         1, 1, 0.8, 0.8
#'     ),
#'     ncol = 4,
#'     byrow = TRUE
#' )
#' ref <- matrix(1, nrow = 20, ncol = 20)
#'
#' reference_use(
#'     target,
#'     samples,
#'     ref,
#'     radius_km = 1000,
#'     k1 = 2,
#'     k2 = 1,
#'     bin_width = 0.1,
#'     interpolate = FALSE,
#'     exclude_slef = FALSE,
#'     num_threads = 1
#' )
reference_use <- function(
        data,
        samples,
        ref_density,
        xy_stats = c(0, 0, 1, 1),
        xy_penalty = 0.0,
        radius_km = 200,
        k1 = 70,
        k2 = 10,
        bin_width = NULL,
        interpolate = TRUE,
        offset = 0,
        confidence = 0.5,
        lambda = 1.0,
        exclude_slef = TRUE,
        drop_features = NULL,
        num_threads = -1,
        weighted_max = FALSE,
        kernel = c("Gaussian", "Cauchy"),
        boost = k2) {

    kernel <- .check_kernel(kernel)
    boost <- .check_boost(boost)
    if (k1 < k2) {
        stop("'k2' must be less than or equal to 'k1'.")
    }
    if (length(confidence) != 1L ||
        !is.finite(confidence) ||
        confidence < 0 ||
        confidence > 1) {
        stop("'confidence' must be one finite number between 0 and 1.")
    }
    if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
        stop("'lambda' must be one finite number greater than zero.")
    }
    if (!is.logical(weighted_max) ||
        length(weighted_max) != 1L ||
        is.na(weighted_max)) {
        stop("'weighted_max' must be one non-missing logical value.")
    }
    if (is.list(samples) && !.is_mat(samples)) {
        stop("Temporal sample lists are not supported by 'reference_use()'.")
    }
    samples <- if (.is_mat(samples)) {
        .check_mat(samples)
    } else {
        stop("'samples' must be a matrix or convertible to one.")
    }

    ref_density <- if (.is_mat(ref_density)) {
        .check_mat(ref_density)
    } else {
        stop("'ref_density' must be a matrix or convertible to one.")
    }
    if (nrow(ref_density) != ncol(ref_density)) {
        warning("Reference density dimensions are not equal!\n")
    }

    if (methods::is(ref_density, "reference_density")) {
        if (is.null(bin_width)) {
            bin_width <- attributes(ref_density)$bin.width
        } else if (bin_width != attributes(ref_density)$bin.width) {
            warning("Provided 'bin_width' differs from reference density attribute.")
        }

        if (is.null(offset)) {
            offset <- attributes(ref_density)$offset
        } else if (offset != attributes(ref_density)$offset) {
            warning("Provided 'offset' differs from reference density attribute.")
        }
    }
    if (is.null(bin_width)) {
        stop("'bin_width' must be supplied when 'ref_density' has no bin-width attribute.")
    }

    if (interpolate) {
        ref_density <- terra::as.matrix(
            terra::disagg(
                terra::rast(.check_mat(unclass(ref_density), name = "ref_density")),
                fact = 2,
                method = "bilinear"
            ),
            wide = TRUE
        )
        bin_width <- bin_width / 2
        offset <- offset * 2
    }
    bin_num <- min(dim(ref_density))

    sample_count <- nrow(samples)
    if (.is_mat(data)) {
        data <- .check_mat(data)
        if (ncol(samples) != ncol(data)) {
            stop("Samples must include all target values and match the columns in 'data'.")
        }

        keep_features <- .keep_rs_features(
            drop_features,
            .num_rs_vars_mat(samples, "samples")
        )
        samples <- .subset_rs_mat(samples, keep_features)
        data <- .subset_rs_mat(data, keep_features)

        result <- reference_use_cpp(
            target_vals = data,
            sample_vals = samples,
            ref_density = ref_density,
            xy_stats = xy_stats,
            xy_penalty = xy_penalty,
            geographic = .is_lonlat(data),
            radius_km = radius_km,
            k_env = k1,
            k_rs = k2,
            bin_width = bin_width,
            bin_num = bin_num,
            offset = offset,
            confidence = confidence,
            boost = boost,
            lambda = lambda,
            exclude_slef = exclude_slef,
            num_threads = num_threads,
            weighted_max = weighted_max,
            kernel = kernel
        )
    } else if (.is_rast(data)) {
        data <- .check_rast(data)
        if (terra::nlyr(data) %% 2L) {
            stop("'data' must contain matching predicted and observed RS layers.")
        }

        if (ncol(samples) == 2L) {
            samples <- cbind(
                samples,
                as.matrix(terra::extract(data, samples, ID = FALSE))
            )
        } else if ((ncol(samples) - 2L) != terra::nlyr(data)) {
            stop("Sample feature count does not match number of raster layers.")
        }

        keep_features <- .keep_rs_features(
            drop_features,
            terra::nlyr(data) / 2L
        )
        samples <- .subset_rs_mat(samples, keep_features)
        data <- .subset_rs_rast(data, keep_features)

        result <- .reference_use_raster(
            data = data,
            samples = samples,
            ref_density = ref_density,
            xy_stats = xy_stats,
            xy_penalty = xy_penalty,
            radius_km = radius_km,
            k1 = k1,
            k2 = k2,
            bin_width = bin_width,
            bin_num = bin_num,
            offset = offset,
            confidence = confidence,
            boost = boost,
            lambda = lambda,
            exclude_slef = exclude_slef,
            num_threads = num_threads,
            weighted_max = weighted_max,
            kernel = kernel
        )
    } else {
        stop("'data' must be a raster, matrix, or convertible object.")
    }

    data.frame(
        id = seq_len(sample_count),
        predicted = result$predicted,
        density = result$density,
        condition = result$condition
    )
}


.reference_use_raster <- function(
        data,
        samples,
        ref_density,
        xy_stats,
        xy_penalty,
        radius_km,
        k1,
        k2,
        bin_width,
        bin_num,
        offset,
        confidence,
        boost,
        lambda,
        exclude_slef,
        num_threads,
        weighted_max,
        kernel) {

    predicted <- numeric(nrow(samples))
    density <- numeric(nrow(samples))
    condition <- numeric(nrow(samples))
    geographic <- .is_lonlat(data)
    block <- terra::blocks(data)

    terra::readStart(data)
    on.exit(terra::readStop(data), add = TRUE)

    for (i in seq_len(block$n)) {
        values <- terra::readValues(
            data,
            row = block$row[i],
            nrows = block$nrows[i],
            mat = TRUE
        )
        first_cell <- (block$row[i] - 1L) * terra::ncol(data) + 1L
        last_cell <- (block$row[i] + block$nrows[i] - 1L) * terra::ncol(data)
        xy <- terra::xyFromCell(data, first_cell:last_cell)

        current <- reference_use_cpp(
            target_vals = cbind(xy, values),
            sample_vals = samples,
            ref_density = ref_density,
            xy_stats = xy_stats,
            xy_penalty = xy_penalty,
            geographic = geographic,
            radius_km = radius_km,
            k_env = k1,
            k_rs = k2,
            bin_width = bin_width,
            bin_num = bin_num,
            offset = offset,
            confidence = confidence,
            boost = boost,
            lambda = lambda,
            exclude_slef = exclude_slef,
            num_threads = num_threads,
            weighted_max = weighted_max,
            kernel = kernel
        )
        predicted <- predicted + current$predicted
        density <- density + current$density
        condition <- condition + current$condition
    }

    list(
        predicted = predicted,
        density = density,
        condition = condition
    )
}
