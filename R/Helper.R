
# are sample points lat-log?
.is_lonlat <- function(x) {
    terra::is.lonlat(
        terra::vect(x),
        perhaps = TRUE,
        warn = FALSE
    )
}


# generate virtual tiles
make_tiles <- function(x, n, rast = FALSE) {
    # x: a path to a raster for creating tiles
    # n: number of tiles
    ## get the number row/col of the tiles
    nc <- ifelse(n %% 2 == 1, 1, floor(sqrt(n / 2)))
    nr <- floor(n / nc)

    r <- terra::rast(x)
    w <- terra::rast(
        terra::ext(r),
        nrows = nr,
        ncols = nc
    )

    if (rast) {
        return(w)
    } else {
        return(
            terra::getTileExtents(r, w)
        )
    }
}


# get raster files
get_rast <- function(x) {
    # x: either raster to character;
    # if only one character supplied, it's considered as a parent path containing the files
    if (is(x, "SpatRaster")) {
        return(x)
    } else if (is(x, "character")) {
        # if length is 1, it means it's a directory path rather than files path
        if (length(x) == 1) {
            r <- terra::rast(
                list.files(x, pattern = "(.tif|.flt)$", full.names = TRUE)
            )
        } else {
            r <- terra::rast(x)
        }
        return(r)
    } else {
        stop("Rasters must either be SpatRaster object or path to the raster layers on disk.")
    }
}

