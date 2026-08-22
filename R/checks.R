# Author: Roozbeh Valavi
# contact: roozbeh.valavi@csiro.au
# Date : Aug-2024
# Version 0.2

# are sample points lat-log?
.is_lonlat <- function(x) {
    return(
        if (inherits(x, "SpatRaster")) {
            terra::is.lonlat(x = x, perhaps = TRUE, warn = FALSE)
        } else if (inherits(x, "matrix")) {
            # Limit to 100 rows to speed up CRS checks
            # Sampling ensures we don’t rely on few points accidentally fall near the projection origin
            rows <- sample(nrow(x), 100, replace = TRUE)
            terra::is.lonlat(x = terra::vect(x[rows, 1:2]), perhaps = TRUE, warn = FALSE)
        } else {
            stop("The 'x' must be a matrix or terra object!")
        }
    )
}


# is matrix or data.frame or data.table
.is_mat <- function(x){
    return(
        inherits(x, c("matrix", "data.table", "data.frame"))
    )
}

# check if it's a matrix if not convert it
.check_mat <- function(x, name = "x") {
    if (.is_mat(x)) {
        if (inherits(x, "matrix")) {
            return(x)
        } else {
            return(
                as.matrix(x)
            )
        }
    } else {
            message(sprintf("'%s' must be a 'matrix', 'data.table', 'data.frame' object.", name))
    }
}

# is it a raster object
.is_rast <- function(x){
    return(
        inherits(
            x,
            c(
                "SpatRaster",
                "RasterStack", "RasterLayer", "RasterBrick",
                "stars",
                "character"
            )
        )
    )
}

# is it a raster or convertible to raster?
.check_rast <- function(r, name = "x"){
    if(!inherits(r, "SpatRaster")){
        tryCatch(
            {
                r <- terra::rast(r)
            },
            error = function(cond) {
                message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
                message(sprintf("'%s' must be a SpatRaster, stars, Raster* object, or path to a raster file on disk.", name))
            }
        )
    }
    return(r)
}

# check for required packages
.check_pkgs <- function(pkg){
    pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
    if(length(pkgna) > 0){
        nm <- paste(pkgna, collapse = ", ")
        message("This function requires these packages: ", nm, "\nWould you like to install them now?\n1: yes\n2: no")
        user <- readline(prompt = paste0("Selection: "))
        if(tolower(user) %in% c("1", "yes", "y")){
            utils::install.packages(pkgna)
        } else{
            stop("Please install these packages for function to work: ", nm)
        }
    }
}


# validate and canonicalise distance-kernel names
.check_kernel <- function(kernel) {
    if (length(kernel) > 1L) {
        kernel <- kernel[[1L]]
    }

    if (!is.character(kernel) ||
        length(kernel) != 1L ||
        is.na(kernel) ||
        !nzchar(kernel)) {
        stop("'kernel' must be 'Gaussian'/'gaussian' or 'Cauchy'/'cauchy'.")
    }

    kernel <- tolower(kernel)
    if (!(kernel %in% c("gaussian", "cauchy"))) {
        stop("'kernel' must be 'Gaussian'/'gaussian' or 'Cauchy'/'cauchy'.")
    }

    kernel
}


# get the number of RS variables from x, y, predicted..., observed... matrix input
.num_rs_vars_mat <- function(x, name = "x") {
    n_vars <- (ncol(x) - 2L) / 2L

    if (ncol(x) < 4L || n_vars != as.integer(n_vars)) {
        stop(
            sprintf(
                "'%s' must contain x, y, predicted RS, and observed RS columns with matching feature counts.",
                name
            )
        )
    }

    as.integer(n_vars)
}


# validate 1-based RS feature indices and return features to keep
.keep_rs_features <- function(drop_features, n_vars) {
    if (!length(drop_features)) {
        return(seq_len(n_vars))
    }

    if (!is.numeric(drop_features) || anyNA(drop_features)) {
        stop("'drop_features' must be an integer vector of RS feature positions.")
    }

    drop_features_int <- as.integer(drop_features)
    if (any(drop_features != drop_features_int)) {
        stop("'drop_features' must contain whole-number RS feature positions.")
    }

    drop_features_int <- sort(unique(drop_features_int))
    if (any(drop_features_int < 1L | drop_features_int > n_vars)) {
        stop(sprintf("'drop_features' must be between 1 and %d.", n_vars))
    }

    setdiff(seq_len(n_vars), drop_features_int)
}


# subset x, y, predicted..., observed... matrix input to the selected RS features
.subset_rs_mat <- function(x, keep_features) {
    n_vars <- .num_rs_vars_mat(x)

    if (length(keep_features) == n_vars) {
        return(x)
    }

    keep_cols <- c(1L, 2L, keep_features + 2L, keep_features + 2L + n_vars)
    x[, keep_cols, drop = FALSE]
}


# subset predicted..., observed... raster layers to the selected RS features
.subset_rs_rast <- function(x, keep_features) {
    n_vars <- terra::nlyr(x) / 2L

    if (length(keep_features) == n_vars) {
        return(x)
    }

    keep_layers <- c(keep_features, keep_features + n_vars)
    x[[keep_layers]]
}


# validate and pack a named list of yearly reference sample matrices
.prepare_temporal_samples <- function(
        samples,
        assessment_year,
        temporal_sigma,
        drop_features = NULL) {

    if (!is.list(samples) || !length(samples)) {
        stop(
            "When 'temporal_sigma' is specified, 'samples' must be a non-empty ",
            "named list of yearly sample matrices."
        )
    }

    sample_year_names <- names(samples)
    if (is.null(sample_year_names) || any(!nzchar(sample_year_names))) {
        stop("Temporal 'samples' must be named with their numeric years.")
    }

    sample_years <- suppressWarnings(as.numeric(sample_year_names))
    if (any(!is.finite(sample_years)) || anyDuplicated(sample_years)) {
        stop("Temporal sample names must be unique numeric years.")
    }
    if (length(assessment_year) != 1L || !is.finite(assessment_year)) {
        stop("'assessment_year' must be one finite numeric value.")
    }
    if (length(temporal_sigma) != 1L ||
        !is.finite(temporal_sigma) ||
        temporal_sigma <= 0) {
        stop("'temporal_sigma' must be one finite number greater than zero.")
    }

    year_order <- order(sample_years)
    sample_years <- sample_years[year_order]
    samples <- samples[year_order]
    samples <- lapply(
        seq_along(samples),
        function(i) {
            if (!.is_mat(samples[[i]])) {
                stop(
                    sprintf(
                        "Temporal sample '%s' must be a matrix or convertible to one.",
                        sample_year_names[year_order][i]
                    )
                )
            }
            .check_mat(samples[[i]])
        }
    )

    n_vars <- .num_rs_vars_mat(samples[[1L]], "samples[[1]]")
    keep_features <- .keep_rs_features(drop_features, n_vars)
    base_xy <- samples[[1L]][, 1:2, drop = FALSE]
    base_pred <- samples[[1L]][, 2L + seq_len(n_vars), drop = FALSE]

    for (i in seq_along(samples)) {
        current <- samples[[i]]
        year_label <- sample_year_names[year_order][i]

        if (.num_rs_vars_mat(current, sprintf("samples[['%s']]", year_label)) != n_vars ||
            nrow(current) != nrow(samples[[1L]])) {
            stop("All temporal sample matrices must have identical dimensions.")
        }
        if (!isTRUE(all.equal(
            current[, 1:2, drop = FALSE],
            base_xy,
            check.attributes = FALSE
        ))) {
            stop("All temporal sample matrices must contain the same XY sites in the same row order.")
        }
        if (!isTRUE(all.equal(
            current[, 2L + seq_len(n_vars), drop = FALSE],
            base_pred,
            check.attributes = FALSE
        ))) {
            stop("Predicted RS values must be constant across temporal sample years.")
        }
        if (anyNA(current)) {
            stop(sprintf("Temporal sample year '%s' contains missing values.", year_label))
        }
    }

    pred_cols <- 2L + keep_features
    obs_cols <- 2L + n_vars + keep_features
    packed <- cbind(
        base_xy,
        samples[[1L]][, pred_cols, drop = FALSE],
        do.call(
            cbind,
            lapply(samples, function(x) x[, obs_cols, drop = FALSE])
        )
    )
    temporal_weights <- exp(
        -0.5 * ((sample_years - assessment_year) / temporal_sigma)^2
    )

    list(
        samples = packed,
        weights = temporal_weights,
        years = sample_years,
        n_vars = n_vars,
        keep_features = keep_features
    )
}

.check_boost <- function(boost) {
    if (is.null(boost) ||
        (length(boost) == 1L && isTRUE(is.na(boost)))) {
        return(NULL)
    }
    if (!is.numeric(boost) ||
        length(boost) != 1L ||
        !is.finite(boost) ||
        boost <= 0) {
        stop("'boost' must be NULL, NA, or one finite number greater than zero.")
    }
    as.numeric(boost)
}
