#' Per-variable departure attribution for HCAS benchmarking
#'
#' Quantifies, for each target location, how much each remote-sensing (RS)
#' variable contributes to the observed departure from its locally selected
#' reference sites, standardised by how tightly those references agree among
#' themselves on that variable. It answers which RS variables drive habitat
#' condition variation, conditional on the local intact ecological distribution,
#' rather than treating importance as a global property of a variable.
#'
#' @details
#' \code{variable_importance()} reuses the same non-temporal three-stage
#' reference selection as \code{\link{benchmark}}: candidates are restricted to
#' those within \code{radius_km}, the \code{k1} nearest in predicted RS
#' space (with an optional XY penalty) are retained, and the \code{k2} with
#' the highest reference-density probability are kept. Each retained reference
#' is weighted with the same selected distance kernel on predicted RS
#' distance used by \code{benchmark()}, so the diagnostic reports on the
#' reference set and weights that benchmarking actually uses. With
#' \code{boost} (the default is \code{k2}), the kernel weight of the retained
#' reference with the highest reference-density probability is multiplied by
#' that factor before signal and noise are computed. Set \code{boost = NULL} or
#' \code{boost = NA} to use the ordinary unboosted kernel weights.
#'
#' For each RS variable \eqn{v} the importance at a target cell is a
#' signal-to-noise ratio
#'
#' \deqn{importance_v = \frac{signal_v}{noise_v + \epsilon}}
#'
#' where, with normalised distance-kernel weights \eqn{w_k} that sum to \eqn{W},
#'
#' \deqn{signal_v = \frac{1}{W} \sum_k w_k\, |obs^{target}_v - obs^{ref,k}_v|}
#' \deqn{noise_v = \frac{1}{W^2} \sum_k \sum_l w_k w_l\, |obs^{ref,k}_v - obs^{ref,l}_v|}
#'
#' The signal is the weighted mean absolute departure of the target from its
#' references on variable \eqn{v}; the noise is the weighted mean absolute
#' difference among the references themselves (a distance-weighted Gini mean
#' difference). Both quantities stay in the L1 geometry of the engine,
#' so the per-variable signals sum to the weighted total observed departure that
#' the reference-density lookup consumes.
#'
#' This is an attribution of the \emph{observed departure} that drives habitat
#' condition, not a sensitivity decomposition of the condition score itself: the
#' departure enters condition non-linearly through the reference-density surface.
#' Importances are univariate and so do not de-correlate variables; collinear RS
#' variables each receive their marginal share. The local-maximum (LDC)
#' component of condition is intentionally excluded, so the metric corresponds to
#' the distance-weighted-mean component of \code{benchmark()}.
#'
#' \code{epsilon} guards the denominator against variables on which the
#' references agree almost perfectly (near-zero noise). Raise it to dampen
#' unstable ratios for low-dispersion variables. Because the noise estimate uses
#' only \code{k2} references, it is itself noisy when the distance weights are
#' concentrated on a single reference; interpret single-reference cells (where
#' noise is zero) with care.
#'
#' Predicted and observed RS variables should be centred and scaled consistently
#' before use, exactly as for \code{\link{benchmark}}, so that per-variable
#' departures and dispersions are comparable across variables.
#'
#' The \code{output} argument selects what each layer/column holds:
#' \describe{
#'   \item{\code{"importance"}}{the signal-to-noise ratio above (default). Best
#'   for ranking and for mapping where a variable is anomalous relative to its
#'   local references. It is open-ended (not normalised) and inherits the
#'   \code{epsilon} sensitivity for low-dispersion variables.}
#'   \item{\code{"signal"}}{the raw weighted absolute departure \eqn{signal_v}.
#'   Per cell these sum to the distance-weighted total observed departure that the
#'   reference-density lookup consumes, so they have a conserved total and are
#'   unaffected by \code{epsilon}.}
#'   \item{\code{"share"}}{the per-cell departure partition
#'   \eqn{signal_v / \sum_u signal_u}. Each cell's variables sum to one (mutually
#'   exclusive), giving "variable \eqn{v} accounts for this fraction of the
#'   observed departure here". Cells whose total departure is zero are
#'   \code{NaN}. This is the most direct percent-attribution map.}
#' }
#'
#' Use \code{\link{aggregate_importance}} to summarise per-cell output across the
#' landscape into a robust ranking (median) and an average attribution (mean of
#' per-cell relative shares).
#'
#' @inheritParams benchmark
#' @param data A matrix, data.frame, or \pkg{terra} \code{SpatRaster} of target
#' RS data. Matrix and data.frame inputs must be organised as \code{x},
#' \code{y}, predicted RS variables, then observed RS variables. Raster inputs
#' must contain predicted RS layers followed by observed RS layers in the same
#' variable order.
#' @param samples A matrix or data.frame of reference sites as \code{x},
#' \code{y}, predicted RS variables, then observed RS variables in the same
#' order as \code{data}. For raster \code{data} it may instead be a two-column
#' coordinate table, in which case raster values are extracted. Temporal sample
#' lists are not supported.
#' @param lambda Positive numeric. Distance-scale bandwidth for the selected
#' \code{kernel} applied to retained-reference predicted RS L1 distances. See
#' \code{\link{benchmark}} for the kernel-specific parameterisation.
#' @param epsilon Numeric, non-negative. Floor added to the noise term to avoid
#' division by (near-)zero reference dispersion. Used only when
#' \code{output = "importance"}.
#' @param output Character. One of \code{"importance"} (default; per-variable
#' signal-to-noise ratio), \code{"signal"} (raw weighted departure
#' contribution), or \code{"share"} (per-cell departure partition summing to one
#' across variables). See Details.
#' @param boost \code{NULL}, \code{NA}, or one positive finite numeric factor.
#' The default \code{k2} multiplies the kernel weight of the
#' highest-probability retained reference by \code{k2} before signal and noise
#' are computed. Use \code{NULL} or \code{NA} for ordinary unboosted kernel
#' weights.
#'
#' @return When \code{data} is a matrix, a numeric matrix with one row per
#' target cell and one column per RS variable. When \code{data} is a raster, a
#' \pkg{terra} \code{SpatRaster} with one layer per RS variable. Cells with
#' missing observed values or no usable references are \code{NaN}; with
#' \code{output = "share"}, cells with zero total departure are also \code{NaN}.
#'
#' @seealso \code{\link{benchmark}}, \code{\link{reference_use}}, and
#' \code{\link{aggregate_importance}}
#' @export
#'
#' @examples
#' target <- matrix(
#'     c(
#'         0, 0, 0.1, 0.2, 0.1, 0.2,
#'         1, 1, 0.8, 0.7, 0.9, 0.6
#'     ),
#'     ncol = 6,
#'     byrow = TRUE
#' )
#' samples <- matrix(
#'     c(
#'         0, 0, 0.1, 0.2, 0.12, 0.18,
#'         1, 1, 0.8, 0.7, 0.85, 0.72,
#'         0, 1, 0.4, 0.5, 0.42, 0.55
#'     ),
#'     ncol = 6,
#'     byrow = TRUE
#' )
#' ref <- matrix(1, nrow = 20, ncol = 20)
#'
#' variable_importance(
#'     target,
#'     samples,
#'     ref,
#'     radius_km = 1000,
#'     k1 = 3,
#'     k2 = 2,
#'     bin_width = 0.1,
#'     interpolate = FALSE,
#'     exclude_slef = FALSE,
#'     num_threads = 1
#' )
variable_importance <- function(
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
        lambda = 1.0,
        epsilon = 1e-6,
        output = c("importance", "signal", "share"),
        exclude_slef = TRUE,
        drop_features = NULL,
        num_threads = -1,
        kernel = c("Gaussian", "Cauchy"),
        boost = k2,
        ...) {

    legacy_k <- intersect(names(list(...)), c("k_pred", "k_obs"))
    if (length(legacy_k)) {
        stop("'k_pred' and 'k_obs' were renamed to 'k1' and 'k2'.")
    }

    output <- match.arg(output)
    kernel <- .check_kernel(kernel)
    boost <- .check_boost(boost)
    # The engine returns either the signal-to-noise importance or the raw
    # per-variable signal; "share" is the per-cell normalised signal computed
    # here, so the engine is asked for "signal" in that case.
    cpp_output <- if (output == "share") "signal" else output

    if (k1 < k2) {
        stop("'k2' must be less than or equal to 'k1'.")
    }
    # Reject benchmark()-only arguments that would otherwise be silently passed
    # through '...' to terra::interpolate and on to the C++ engine, which does
    # not accept them. 'confidence' is the most common copy-paste mistake.
    benchmark_only <- intersect(
        names(list(...)),
        c("confidence", "make_su", "temporal_correct",
          "assessment_year", "temporal_sigma")
    )
    if (length(benchmark_only)) {
        stop(
            "variable_importance() does not accept benchmark() argument(s): ",
            paste(benchmark_only, collapse = ", "), ".\n",
            "  In particular, 'confidence' has no effect here: the metric ",
            "attributes the distance-weighted-mean departure and deliberately ",
            "excludes the LDC/confidence component. Remove it from the call."
        )
    }
    if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
        stop("'lambda' must be one finite number greater than zero.")
    }
    if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0) {
        stop("'epsilon' must be one finite, non-negative number.")
    }
    if (is.list(samples) && !.is_mat(samples)) {
        stop("Temporal sample lists are not supported by 'variable_importance()'.")
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
        var_names <- .rs_var_names(data, keep_features)

        result <- variable_importance_cpp(
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
            lambda = lambda,
            epsilon = epsilon,
            output = cpp_output,
            exclude_slef = exclude_slef,
            num_threads = num_threads,
            kernel = kernel,
            boost = boost
        )
        if (output == "share") {
            result <- .importance_shares(result)
        }
        colnames(result) <- var_names
        return(result)

    } else if (.is_rast(data)) {
        .check_pkgs("terra")
        data <- .check_rast(data)
        if (terra::nlyr(data) %% 2L) {
            stop("'data' must contain matching predicted and observed RS layers.")
        }

        if (ncol(samples) == 2L) {
            cat("Extracting sample values...\n")
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
        var_names <- .rs_var_names_rast(data)

        out_rast <- terra::interpolate(
            object = data,
            model = list(),
            fun = .variable_importance_predict,
            var_names = var_names,
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
            lambda = lambda,
            epsilon = epsilon,
            output = cpp_output,
            exclude_slef = exclude_slef,
            num_threads = num_threads,
            kernel = kernel,
            boost = boost,
            ...
        )
        if (output == "share") {
            total <- terra::app(out_rast, "sum")
            total <- terra::ifel(total == 0, NA, total)
            out_rast <- out_rast / total
        }
        names(out_rast) <- var_names
        return(out_rast)

    } else {
        stop("'data' must be a raster, matrix, or convertible object.")
    }
}


# terra::interpolate prediction helper: returns one column per RS variable.
# Surfaces the real error rather than silently NA-ing the block: the C++ engine
# returns NaN rows for empty/missing cells itself, so any thrown error is a
# structural problem that would affect every block and must abort the run.
.variable_importance_predict <- function(model, newdata, var_names, ...) {
    dat <- as.matrix(newdata)

    out <- tryCatch(
        variable_importance_cpp(target_vals = dat, ...),
        error = function(cond) {
            stop(
                "variable_importance() failed while evaluating a raster block: ",
                conditionMessage(cond),
                call. = FALSE
            )
        }
    )

    colnames(out) <- var_names
    out
}


# observed RS variable names from a x,y,pred...,obs... matrix
.rs_var_names <- function(x, keep_features) {
    n_vars <- (ncol(x) - 2L) / 2L
    obs_names <- colnames(x)[(2L + n_vars + 1L):ncol(x)]
    if (is.null(obs_names) || any(!nzchar(obs_names))) {
        return(paste0("var", keep_features))
    }
    obs_names
}


# observed RS variable names from a pred...,obs... raster
.rs_var_names_rast <- function(x) {
    n_vars <- terra::nlyr(x) / 2L
    obs_names <- names(x)[(n_vars + 1L):terra::nlyr(x)]
    if (is.null(obs_names) || any(!nzchar(obs_names))) {
        return(paste0("var", seq_len(n_vars)))
    }
    obs_names
}


# per-cell normalisation of a per-variable signal matrix into shares summing to
# one; rows whose total signal is zero or non-finite become NA
.importance_shares <- function(x) {
    rs <- rowSums(x)
    rs[!is.finite(rs) | rs == 0] <- NA_real_
    x / rs
}


#' Summarise per-variable importance across the landscape
#'
#' Aggregates the per-cell output of \code{\link{variable_importance}} into a
#' landscape-level characterisation of which RS variables most drive habitat
#' condition variation.
#'
#' @details
#' Two complementary summaries are returned for each variable. The
#' \code{median} importance gives a robust ranking that is insensitive to a
#' minority of cells with extreme ratios. The \code{mean_share} is the mean over
#' cells of the variable's relative share of importance,
#' \eqn{importance_v / \sum_u importance_u}, expressed as a proportion; multiply
#' by 100 for an average percent attribution. The two can disagree when a
#' variable is usually modest but occasionally dominant, so reporting both is
#' recommended.
#'
#' Shares are computed only over cells where every variable's importance is
#' finite and the row sum is positive. Because importances are conditional on
#' overlapping local reference sets, nearby cells are not independent draws;
#' treat the aggregate as a description of the landscape rather than a basis for
#' naive inferential standard errors. Where the premise that importance is local
#' matters, aggregate within ecologically meaningful strata (for example region
#' or vegetation class) instead of a single global summary.
#'
#' @param x A numeric matrix (one row per cell, one column per variable) as
#' returned by \code{\link{variable_importance}} for matrix input, or a
#' \pkg{terra} \code{SpatRaster} as returned for raster input. Raster values are
#' read into memory, so summarise per tile for very large analyses.
#'
#' @return A data.frame with one row per variable, ordered by descending
#' \code{median}, containing:
#' \itemize{
#'   \item \code{variable}: variable name.
#'   \item \code{median}: median per-cell importance (robust ranking).
#'   \item \code{mean_share}: mean per-cell relative share (average attribution).
#'   \item \code{rank}: 1-based rank by \code{median}.
#'   \item \code{n}: number of cells contributing to \code{mean_share}.
#' }
#'
#' @seealso \code{\link{variable_importance}}
#' @export
#'
#' @examples
#' imp <- matrix(
#'     c(
#'         0.8, 0.2, 0.1,
#'         0.6, 0.3, 0.1,
#'         0.7, 0.2, 0.2
#'     ),
#'     ncol = 3,
#'     byrow = TRUE,
#'     dimnames = list(NULL, c("a", "b", "c"))
#' )
#' aggregate_importance(imp)
aggregate_importance <- function(x) {
    if (inherits(x, "SpatRaster")) {
        var_names <- names(x)
        x <- terra::values(x, mat = TRUE)
    } else if (.is_mat(x)) {
        x <- .check_mat(x)
        var_names <- colnames(x)
    } else {
        stop("'x' must be a matrix or SpatRaster of per-cell importances.")
    }

    if (is.null(var_names) || any(!nzchar(var_names))) {
        var_names <- paste0("var", seq_len(ncol(x)))
    }

    med <- apply(x, 2L, stats::median, na.rm = TRUE)

    row_sums <- rowSums(x)
    valid <- is.finite(row_sums) & row_sums > 0
    if (any(valid)) {
        shares <- x[valid, , drop = FALSE] / row_sums[valid]
        mean_share <- colMeans(shares, na.rm = TRUE)
        n_valid <- sum(valid)
    } else {
        mean_share <- rep(NA_real_, ncol(x))
        n_valid <- 0L
    }

    out <- data.frame(
        variable = var_names,
        median = as.numeric(med),
        mean_share = as.numeric(mean_share),
        n = n_valid,
        stringsAsFactors = FALSE
    )
    out <- out[order(-out$median), , drop = FALSE]
    out$rank <- seq_len(nrow(out))
    rownames(out) <- NULL
    out[, c("variable", "median", "mean_share", "rank", "n")]
}
