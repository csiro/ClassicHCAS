pkgname <- "ClassicHCAS"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ClassicHCAS')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("benchmark")
### * benchmark

flush(stderr()); flush(stdout())

### Name: benchmark
### Title: Condition benchmarking of target points
### Aliases: benchmark

### ** Examples




cleanEx()
nameEx("calibrate")
### * calibrate

flush(stderr()); flush(stdout())

### Name: calibrate
### Title: Calibrate habitat condition output
### Aliases: calibrate

### ** Examples




cleanEx()
nameEx("normalise")
### * normalise

flush(stderr()); flush(stdout())

### Name: normalise
### Title: Clean and normalise HCAS reference density
### Aliases: normalise

### ** Examples




cleanEx()
nameEx("radial_count")
### * radial_count

flush(stderr()); flush(stdout())

### Name: radial_count
### Title: Number of samples within a radius
### Aliases: radial_count

### ** Examples




cleanEx()
nameEx("ref_density")
### * ref_density

flush(stderr()); flush(stdout())

### Name: ref_density
### Title: Reference density surface
### Aliases: ref_density

### ** Examples




cleanEx()
nameEx("tiling")
### * tiling

flush(stderr()); flush(stdout())

### Name: tiling
### Title: Generate raster tiles using raster or matrix data
### Aliases: tiling

### ** Examples

## Not run: 
##D library(ClassicHCAS)
##D library(terra)
##D 
##D r <- rast(nrows=100, ncols=100)
##D values(r) <- runif(ncell(r))
##D 
##D # Balanced tiles
##D tiles_poly <- ClassicHCAS::tiling(r, n_tiles = 4, balanced = TRUE, spatial = TRUE)
##D 
##D # Rectangular tiles
##D tiles_rect <- ClassicHCAS::tiling(r, n_tiles = 4, balanced = FALSE)
##D 
##D # Balanced tiles from a matrix
##D mat <- as.matrix(r)
##D ext <- ext(r)
##D tiles_from_matrix <- ClassicHCAS::tiling(mat, n_tiles = 4, balanced = TRUE, extent = ext)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
