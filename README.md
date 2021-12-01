___
output: github_content
___



# fastMCD

## Introduction
The minimum covariance determinant (MCD) method obtains robust estimates of location and scale for multivariate data by identifying a subset of _h_ observations with the lowest determinant in its covariance matrix. This method is highly robust and accurate but is infeasible for large datasets. The FAST-MCD algorithm, described by Rousseeuw and van Driessen (1999), approximates the MCD estimate with an iterative procedure applied to a set of initial subsamples. FAST-MCD is far faster than MCD and typically is accurate even for large datasets. Superior implementations exist in R already, however the `fastMCD` package represents my attempt to independently implement FAST-MCD and understand the underlying algorithm.


## Installation

```r
devtools::install_github("frankp-0/fastMCD")
## Skipping install of 'fastMCD' from a github remote, the SHA1 (ebb2af9b) has not changed since last install.
##   Use `force = TRUE` to force installation
```

## Usage

```r
library(fastMCD)
data(smallX)
```
