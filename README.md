# ClassicHCAS

[![R build
status](https://github.com/csiro/ClassicHCAS/workflows/R-CMD-check/badge.svg)](https://github.com/csiro/ClassicHCAS/actions)
![maintenance-status](https://img.shields.io/badge/maintenance-active-brightgreen.svg)

## What is ClassicHCAS?

`ClassicHCAS` is the R implementation of the Habitat Condition Assessment System (HCAS), a CSIRO method for estimating habitat condition from remote sensing data and reference ecosystem samples. HCAS compares the relationship between predicted and observed remote sensing variables at reference sites with the same relationship at target locations. Sites that behave more like the reference system receive higher habitat condition scores.

At a high level, HCAS is a reference-based habitat condition estimation method. It uses relatively intact sites to learn the expected relationship between modelled and observed remote sensing signals, then measures how closely other locations follow that reference pattern. The result is a spatially explicit condition surface that can be used in habitat assessment, monitoring, and large-area reporting.

This package consolidates the original HCAS toolset into a single R package and uses `Rcpp` for the computationally intensive steps.

## What the package does

`ClassicHCAS` supports the main HCAS workflow:

- `ref_density()` builds a reference density surface from reference-site data and remote sensing information.
- `normalise()` trims and normalises that surface for benchmarking.
- `benchmark()` scores target sites or rasters against the reference relationship.
- `calibrate()` rescales raw condition values to a final 0 to 1 condition scale.

It also includes operational helpers:

- `radial_count()` counts nearby reference samples for each raster cell.
- `tiling()` creates balanced or rectangular tiles for large processing runs.

## Expected inputs

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

## Typical workflow

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
  x_values = c(0, 0.001, 0.03, 0.038),
  y_values = c(0, 0.2, 0.8, 1)
)
```

## Australian application and citation

For the Australian national application of HCAS, including the official project overview, product guides, factsheets, and current data releases, see the CSIRO HCAS project page:

- <https://research.csiro.au/biodiversity-knowledge/projects/hcas/>

That page is specifically about the Australian national application of HCAS. CSIRO describes it as Australia's first consistent, repeatable, and cost-efficient national habitat condition assessment, monitoring, and reporting capability built on Earth observation products. It is intended to support management prioritisation, national environmental reporting, and analysis of natural and non-natural influences on habitat condition across continental Australia.

If you are citing the current Australian data collection referenced on that page, use the HCAS 3.3 citation:

- Valavi R, Levick SR, Lehmann EA, Liu N, Giljohann KM, Williams KJ, Collings S, Johnson S, Botha EJ, Munroe SEM, Van Niel TG, Newnham G, Paget M, Malley C, Carlile P, Gunawardana D, Lyon P, Richards AE, Tetreault Campbell S and Ferrier S (2025c). *HCAS 3.3 (1988-2024) base model estimate of habitat condition (90m grid), National Connectivity Index 2.0 (NCI) and annual time series for continental Australia*. Data collection 65549. CSIRO, Canberra, Australia. Citation text from the CSIRO HCAS project page: <https://research.csiro.au/biodiversity-knowledge/projects/hcas/>

For the latest Australian release status, product guides, and download links, use the CSIRO project page above as the canonical source.
