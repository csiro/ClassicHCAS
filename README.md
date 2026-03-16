# ClassicHCAS

[![R build
status](https://github.com/csiro/ClassicHCAS/workflows/R-CMD-check/badge.svg)](https://github.com/csiro/ClassicHCAS/actions)
![maintenance-status](https://img.shields.io/badge/maintenance-active-brightgreen.svg)

## What Is ClassicHCAS?

`ClassicHCAS` is the R implementation of the Habitat Condition Assessment System (HCAS), a CSIRO method for estimating habitat condition from remote sensing data and reference ecosystem samples.

HCAS compares the relationship between predicted and observed remote sensing variables at reference sites with the same relationship at target locations. Sites that behave more like the reference system receive higher habitat condition scores.

This package consolidates the original HCAS toolset into a single R package and uses `Rcpp` for the computationally intensive steps.

## What The Package Does

`ClassicHCAS` supports the main HCAS workflow:

- `ref_density()` builds a reference density surface from reference-site data.
- `normalise()` trims and normalises that surface for benchmarking.
- `benchmark()` scores target sites or rasters against the reference relationship.
- `calibrate()` rescales raw condition values to a final 0 to 1 condition scale.

It also includes operational helpers:

- `radial_count()` counts nearby reference samples for each raster cell.
- `tiling()` creates balanced or rectangular tiles for large processing runs.

## Expected Inputs

Most functions accept either:

- matrices or data frames arranged as `x`, `y`, predicted remote sensing variables, and observed remote sensing variables
- `terra::SpatRaster` objects with predicted and observed layers in the same order

Predicted and observed variables should be centred and scaled consistently before creating the reference density or benchmarking target sites.

## Installation

`ClassicHCAS` is not currently on CRAN. Install it from GitHub:

```r
install.packages("remotes")
remotes::install_github("csiro/ClassicHCAS")
```

Building from source requires a working toolchain. On Windows this usually means [Rtools](https://cran.r-project.org/bin/windows/Rtools/). On macOS, OpenMP support may be needed for faster parallel execution.

## Typical Workflow

```r
library(ClassicHCAS)

ref <- ref_density(
  data = reference_data,
  radius_km = 1000,
  bin_width = 0.05
)

ref_norm <- normalise(ref)

condition <- benchmark(
  data = target_data,
  samples = benchmark_samples,
  ref_density = ref_norm
)

condition_calibrated <- calibrate(
  condition,
  x_values = c(0, 0.2, 0.7, 1)
)
```

## What HCAS Means In Practice

At a high level, HCAS is a reference-based habitat condition method. It uses relatively intact sites to learn the expected relationship between modelled and observed remote sensing signals, then measures how closely other locations follow that reference pattern. The result is a spatially explicit condition surface that can be used in habitat assessment, monitoring, and large-area reporting.
