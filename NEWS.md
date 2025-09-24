# Version 1.0.0
* Complete rewrite of the main functions (`histogram` and `benchmark`) using vectorised operations, resulting in more than a 10× speed improvement.
* The distance calculations are now changed from double to float32, providing significant speed improvements with no loss of accuracy.
* The histogram now performs only one-way pairwise distance calculations. As a result, the raw histogram values are exactly halved compared to before. This has no impact on the normalised histogram or the final output.
* The "corner value” in histogram (previously calculated as the count of values) is no longer computed, since it is simply equal to the number of samples.
* The `calibrate` function no longer performs interpolation. The output is now fully fitted using a monotonic spline.
* The NAN pixels are now directly handled in the C++ side.

# Version 0.2.0
* Added `drop_features` parameter to fully exclude specific remote sensing variables from both the `histogram` and `benchmark` functions.
* Added a condition to ensure that `k_obs`  is less than or equal to `k_pred`.
* A new histogram normalisation method has been implemented in R using the `legacy = FALSE` argument to mitigate edge effects.
* The zero-zero point (reference self-count) in the raw histogram is now excluded prior to normalisation.
* The `histogram` function arguments have changed to `data` and `samples`.

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
