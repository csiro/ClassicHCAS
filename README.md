# ClassicHCAS
## The Classic Habitat Condition Assessment System (HCAS)

The `ClassicHCAS` package provides a comprehensive suite of functions for assessing habitat condition using remote sensing indices and reference ecosystem samples. This package is a complete re-write of CSIRO’s earlier software collection, **HCAS**, originally developed by Tom Harwood. All the tools and functions from the previous HCAS software are now available in this R package, which supports a variety of data formats. 
The calculation of HCAS products across the entire continent is now fully automated through the [hcas-workflow](https://bitbucket.csiro.au/projects/HCAS3/repos/hcas_workflow/browse) Python repository, using the `ClassicHCAS` package.

**Note**: The `ClassicHCAS` package is still under active development, so the user interface might completely change.

## Installation

The `Rtools` toolchain is required to install/build the package on Windows, see [the following link](https://cran.r-project.org/bin/windows/Rtools/) for downloading and installing `Rtools`.

To install the package on a CSIRO system use:

```r
remotes::install_git("https://bitbucket.csiro.au/scm/hcas3/classichcas.git")
```

This might ask your Bitbcket username and password.

## Usage

Examples will be provided here soon.
