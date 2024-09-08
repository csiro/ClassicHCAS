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
