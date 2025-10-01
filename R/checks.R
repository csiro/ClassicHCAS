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
