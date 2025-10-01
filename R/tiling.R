#' Generate raster tiles using raster or matrix data
#'
#' Creates either rectangular or balanced tiles for raster or matrix data. Balanced tiles
#' attempt to distribute the data values evenly across tiles, while rectangular tiles simply
#' divide the raster into equal-sized grids. This function helps creating tiles for running
#' \code{benchmark} function over multiple tiles or systems (e.g. in a cluster). The output
#' of the \code{\link{proximity}} function can be used to balance the run time over each tile.
#'
#' @param data A `SpatRaster` object or a numeric matrix. For balanced tiling, either a matrix
#' or raster can be used. For rectangular tiling, a raster is required.
#' @param n_tiles Integer. Number of tiles to generate.
#' @param balanced Logical. If `TRUE` (default), tiles are created to balance the sum of values within each tile.
#'   If `FALSE`, rectangular tiles of equal size are generated.
#' @param method Character. One of `"best"`, `"row"`, `"col"`, or `"both"` (default `"best"`).
#'   Specifies how the balanced tiles should be split:
#'   \describe{
#'     \item{"best"}{Automatically chooses the split direction that minimizes imbalance.}
#'     \item{"row"}{Always splits by rows.}
#'     \item{"col"}{Always splits by columns.}
#'     \item{"both"}{Splits along both dimensions, favoring the longer dimension.}
#'   }
#' @param exact Logical. If `TRUE` (default), ensures exactly `n_tiles` are produced in balanced mode.
#' @param weighted Logical. If `TRUE` (default), tile weights are based on the scaled data values.
#'   If `FALSE`, all non-zero values are treated equally.
#' @param spatial Logical. If `TRUE`, returns a `SpatVector` polygon of tiles.
#'   If `FALSE` (default), returns a matrix or data frame with tile extents.
#' @param extent Optional. A `terra::ext` object specifying the raster extent. Required when `data` is a matrix.
#'
#' @return Either:
#' \itemize{
#'   \item A `SpatVector` of polygons representing tiles (if `spatial = TRUE`), or
#'   \item A matrix/data.frame with columns `"xmin"`, `"xmax"`, `"ymin"`, `"ymax"` for each tile (if `spatial = FALSE`).
#' }
#'
#' @details
#' - Rectangular tiling divides the raster into equal-sized tiles regardless of data values.
#' - Balanced tiling attempts to split the data such that the sum of the cell values in each tile is roughly equal.
#' - The splitting process recursively divides the matrix/raster along rows or columns, according to the `method`.
#' - NA values are treated as zero for the purpose of balanced tiling.
#'
#' @examples
#' \dontrun{
#' library(ClassicHCAS)
#' library(terra)
#'
#' r <- rast(nrows=100, ncols=100)
#' values(r) <- runif(ncell(r))
#'
#' # Balanced tiles
#' tiles_poly <- ClassicHCAS::tiling(r, n_tiles = 4, balanced = TRUE, spatial = TRUE)
#'
#' # Rectangular tiles
#' tiles_rect <- ClassicHCAS::tiling(r, n_tiles = 4, balanced = FALSE)
#'
#' # Balanced tiles from a matrix
#' mat <- as.matrix(r)
#' ext <- ext(r)
#' tiles_from_matrix <- ClassicHCAS::tiling(mat, n_tiles = 4, balanced = TRUE, extent = ext)
#' }
#'
#' @seealso
#' \code{\link[terra]{rast}}, \code{\link[terra]{as.polygons}}, \code{\link[terra]{ext}}
#'
#' @export
tiling <- function(
        data,            # matrix (only for balanced) or raster for creating the tiles
        n_tiles,         # number of tiles
        balanced = TRUE, # should it be balanced tiles, otherwise rectangular
        method = c("best", "row", "col", "both"), # balanced spliting method
        exact = TRUE,    # the exact number of tiles for balanced tiles?
        weighted = TRUE, # weights must be non-zero or NA; the NA have no weights for tiling
        spatial = FALSE, # return a polygon, otherwise a csv
        extent = NULL    # enforce raster extent to make sure it won't miss any pixels
) {
    # define the choice of splitting
    method <- match.arg(method[1], choices = c("best", "row", "col", "both"))

    equi_tiles <- function(r, n, sp = FALSE) {
        nc <- if (n %% 2 == 1) 1 else floor(sqrt(n / 2))
        nr <- floor(n / nc)

        w <- terra::rast(terra::ext(r), nrows = nr, ncols = nc)

        if (sp) {
            return(terra::as.polygons(w))
        } else {
            return(terra::getTileExtents(r, w))
        }
    }

    if (!balanced) {
        if (any(methods::is(data, "SpatRaster"))) {
            return(
                equi_tiles(
                    r = data,
                    n = n_tiles,
                    sp = spatial
                )
            )
        } else {
            stop("Equal-siezed tiles only works with rasters.")
        }
    }

    # get data and make sure there's no NA
    if (any(methods::is(data, "SpatRaster"))) {
        x <- terra::as.matrix(data, wide = TRUE)
    } else if (any(is(data, "matrix"))) {
        x <- data
    } else {
        stop("'data' must be a matrix or SpatRaster!")
    }

    # make sure every pixel is counted
    x <- ifelse(is.na(x), 0, scales::rescale(x, to = c(1, 100)))

    # make all cells weighted equally
    if (!weighted) {
        x <- ifelse(x > 0, 1, 0)
    }

    # find the best splitting point
    split_point <- function(x, byrow = TRUE, diff = FALSE) {
        dimsum <- if (byrow) rowSums(x, na.rm = TRUE) else colSums(x, na.rm = TRUE)
        ave <- sum(dimsum, na.rm = TRUE) / 2
        cs <- cumsum(dimsum) - ave
        pt <- which.min(abs(cs))
        # just return the difference in points?
        if (diff) {
            return(
                abs(
                    sum(dimsum[1:pt], na.rm = TRUE) - sum(dimsum[(pt+1):length(dimsum)], na.rm = TRUE)
                )
            )
        }

        return(
            pt
        )
    }
    # split by row or column
    by_rows <- function(x, by = "best") {
        if (by == "best") {
            brow <- split_point(x = x, byrow = TRUE, diff = TRUE)
            bcol <- split_point(x = x, byrow = FALSE, diff = TRUE)
            # if difference in node stats are higher in by-column, then choose by-rows
            return(
                ifelse(bcol >= brow, TRUE, FALSE)
            )
        } else if (by == "row") {
            return(TRUE)
        } else if (by == "col") {
            return(FALSE)
        } else {
            return(
                ifelse(nrow(x) >= ncol(x), TRUE, FALSE)
            )
        }
    }
    # filter the matrix based on lookup
    filter_mat <- function(mat, lookup, i) {
        r1 <- lookup[i, "r1"]
        r2 <- lookup[i, "r2"]
        c1 <- lookup[i, "c1"]
        c2 <- lookup[i, "c2"]
        return(
            mat[r1:r2, c1:c2, drop = FALSE]
        )
    }
    # calculate the sums
    sum_mat <- function(mat, lookup, i) {
        r1 <- lookup[i, "r1"]
        r2 <- lookup[i, "r2"]
        c1 <- lookup[i, "c1"]
        c2 <- lookup[i, "c2"]
        return(
            sum(mat[r1:r2, c1:c2], na.rm = TRUE)
        )
    }

    n_rows <- nrow(x)
    n_cols <- ncol(x)
    result <- matrix(1, n_rows, n_cols)

    total_sum <- sum(x, na.rm = TRUE)
    target_sum <- total_sum / n_tiles

    # initial row/col id
    row_start <- 1
    col_start <- 1
    row_end <- n_rows
    col_end <- n_cols
    id <- 1

    lookup <- data.frame(
        id = 1,
        r1 = row_start,
        r2 = row_end,
        c1 = col_start,
        c2 = col_end,
        wt = total_sum
    )

    # find the splits
    while (max(lookup$wt) > target_sum && ifelse(exact, nrow(lookup) < n_tiles, TRUE)) {
        r <- which.max(lookup$wt)
        sub_mat <- filter_mat(x, lookup, r)
        split_row <- by_rows(sub_mat, by = method)
        # update look up table
        if (split_row) {
            id <- id + 1
            pnt <- split_point(x = sub_mat, byrow = TRUE)
            new_end_row <- (lookup[r, "r1"] + pnt - 1)
            row_start <- min(new_end_row + 1, n_rows)
            row_end <- lookup[r, "r2"] # first update this
            col_start <- lookup[r, "c1"]
            col_end <- lookup[r, "c2"]
            lookup[r, "r2"] <- new_end_row
            lookup[nrow(lookup) + 1, ] <- c(id, row_start, row_end, col_start, col_end, 0)
            lookup[nrow(lookup), "wt"] <- sum_mat(x, lookup, nrow(lookup))
            lookup[r, "wt"] <- sum_mat(x, lookup, r)
        } else {
            id <- id + 1
            pnt <- split_point(x = sub_mat, byrow = FALSE)
            new_end_col <- (lookup[r, "c1"] + pnt - 1)
            row_start <- lookup[r, "r1"]
            row_end <- lookup[r, "r2"]
            col_start <- min(new_end_col + 1, n_cols)
            col_end <- lookup[r, "c2"]
            lookup[r, "c2"] <- new_end_col
            lookup[nrow(lookup) + 1, ] <- c(id, row_start, row_end, col_start, col_end, 0)
            lookup[nrow(lookup), "wt"] <- sum_mat(x, lookup, nrow(lookup))
            lookup[r, "wt"] <- sum_mat(x, lookup, r)
        }

        result[row_start:row_end, col_start:col_end] <- id
    }

    # get the ext and make it a polygon
    if (any(methods::is(data, "SpatRaster"))) {
        if (is.null(extent)) {
            extent <- terra::ext(data)
        }
    } else {
        if (is.null(extent)) {
            stop("For matrix 'data' you need to provide the 'extent'.")
        }
    }

    outpoly <- terra::as.polygons(
        terra::rast(result, extent = extent)
    )


    if (spatial) {
        if (any(methods::is(data, "SpatRaster"))) {
            terra::crs(outpoly) <- terra::crs(data)
        }

        return(
            outpoly
        )
    } else {
        return(
            .get_tile_extent(outpoly)
        )
    }
}

# get extent of each tile
.get_tile_extent <- function(x) {
    nr <- nrow(x)
    out <- matrix(0, nrow = nr, ncol = 4)
    colnames(out) <- c("xmin", "xmax", "ymin", "ymax")
    for(i in 1:nr) {
        out[i, ] <- terra::ext(x[i, ])[1:4]
    }

    return(
        out
    )
}

