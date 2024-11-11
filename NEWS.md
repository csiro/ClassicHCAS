# Version 0.1.6
* Added a `scale_factor` parameter to allow user-defined correction of geographic CRS distance calculations, enhancing flexibility in handling distance conversions (previously set to a default value only).
* Enhanced C++ code for improved efficiency in point class creation.

# Version 0.1.5
* Replaced the natural spline in the `calibrate` function with a monotonically increasing spline function.

# Version 0.1.4
* Deprecated the linear (piece-wise) calibration method, superseded with the spline method.
* The arguments of `calibrate` are updated with `x_values` and `y_values`.
* Fixed the benchmarking `self_exclude` floating point error.

# Version 0.1.3
* Added the `proximity` function to count the number of samples within a specified radius.
* Removed the requirement for `add_xy` in the `benchmark()` function.
* Removed `filename` and `wopt` parameters from all functions and replaced them with additional arguments `...` for more flexibility.
* Internal `terra` functions have been improved by eliminating the need for explicit creation of x and y coordinates.

# Version 0.1.2
* Check for availability of the `terra` package is added when input is a raster
* improved documentation and imports
* a bug fix in benchmarking

# Version 0.1.1
* Added the geographic distance penalty for selecting the nearest predicted neighbours (`xy_penalty` and `xy_stats` parameters)
* Added an option for excluding self-assessment for benchmark sample sites (default).

# Version 0.1.0
* Translated all legacy HCAS codes to R and Rcpp
