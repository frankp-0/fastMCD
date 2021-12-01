
# fastMCD

## Introduction

The minimum covariance determinant (MCD) method obtains robust estimates
of location and scale for multivariate data by identifying a subset of
*h* observations with the lowest determinant in its covariance matrix.
This method is highly robust and accurate but is infeasible for large
datasets. The FAST-MCD algorithm, described by Rousseeuw and van
Driessen (1999), approximates the MCD estimate with an iterative
procedure applied to a set of initial subsamples. FAST-MCD is far faster
than MCD and typically is accurate even for large datasets. Superior
implementations exist in R already, however the `fastMCD` package
represents my attempt to independently implement FAST-MCD and understand
the underlying algorithm.

## Installation

``` r
devtools::install_github("frankp-0/fastMCD")
## Downloading GitHub repo frankp-0/fastMCD@HEAD
##      checking for file ‘/private/var/folders/sd/f9yrf8zj0vbdv0ny3rzn00pm0000gn/T/RtmpPsb7uO/remotes70c07cf6325b/frankp-0-fastMCD-e9054c3/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/sd/f9yrf8zj0vbdv0ny3rzn00pm0000gn/T/RtmpPsb7uO/remotes70c07cf6325b/frankp-0-fastMCD-e9054c3/DESCRIPTION’
##   ─  preparing ‘fastMCD’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
## ─  looking to see if a ‘data/datalist’ file should be added
##        NB: this package now depends on R (>= 3.5.0)
##        WARNING: Added dependency on R >= 3.5.0 because serialized objects in
##      serialize/load version 3 cannot be read in older versions of R.
##      File(s) containing such objects:
##        ‘fastMCD/data/bigX.rda’ ‘fastMCD/data/smallX.rda’
##   ─  building ‘fastMCD_0.0.0.9000.tar.gz’
##      
## 
```

## Usage

``` r
library(fastMCD)
data(smallX)
```
