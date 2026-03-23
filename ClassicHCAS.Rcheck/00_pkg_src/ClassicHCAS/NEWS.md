# Version 1.1.0
* Consolidated package colour helpers into a single `palettes()` function and removed `hcas_color()` and `ref_density_color()`.
* Changed `drop_features` so excluded RS variables are removed from the active feature set before C++ distance calculations, instead of being zeroed in place.

# Version 1.0.0
* Renamed the reference density function to `ref_density()`, renamed the radial counting function (from `proximity()`) to `radial_count()`, and changed the reference density object class to `reference_density`.
* The order of input data has completely changed and both `ref_density()` and `benchmark()` functions require data with **x**, **y**, **predicted**, and **observed** remote sensing variable order.
* Complete re-write of the main C++ functions (for `ref_density()`, `benchmark()`, and `radial_count()`) using vectorised operations, including changing feature space distance calculations from `double` to `float32`, resulting in significant speed improvement (~10x) with no loss of accuracy.
* The Eigen C++ library is adapted as the main matrix engine.
* The dependency on KDtrees is dropped while keeping or even improving the speed.
* The reference density now performs only one-way pairwise distance calculations. As a result, the raw reference density values are exactly halved compared to before. This has no impact on the normalised reference density or the final output.
* The "corner value" in the reference density (previously calculated as the count of values) is no longer computed, since it is simply equal to the number of samples.
* The `calibrate()` function no longer performs interpolation. The output is now fully fitted using a monotonic spline.
* A fast spatial distance calculation is implemented (no difference for projected CRS with 0.01 meters accuracy).
* The geographic distance for lat/long is now corrected for the latitude of the source cell, not just the radius transformation to degrees.
* A new function (`tiling()`) is added for making raster tiles using raster or matrix data.
* The `NaN` pixels are now directly handled within C++ code.
* Expanded the test suite with additional C++-focused coverage, including edge cases and thread-consistency checks for core workflows.

# Version 0.2.0
* Added `drop_features` parameter to fully exclude specific remote sensing variables from both the `ref_density` and `benchmark` functions.
* Added a condition to ensure that `k_obs`  is less than or equal to `k_pred`.
* A new reference density normalisation method has been implemented in R using the `legacy = FALSE` argument to mitigate edge effects.
* The zero-zero point (reference self-count) in the raw reference density is now excluded prior to normalisation.
* The `ref_density` function arguments have changed to `data` and `samples`.

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
