---
output: github_document
---



# fastMCD

## Overview
`fastMCD` performs robust mean and covariance estimation using the FAST-MCD algorithm described by Rousseeuw and van Driessen (1999).

## Introduction
The minimum covariance determinant (MCD) method obtains robust estimates of location and scale for multivariate data by identifying a subset of _h_ observations with the lowest determinant in its covariance matrix. This method is highly robust and accurate but is infeasible for large datasets. The FAST-MCD algorithm, described by Rousseeuw and van Driessen (1999), approximates the MCD estimate with an iterative procedure applied to a set of initial subsamples. FAST-MCD is far faster than MCD and typically is accurate even for large datasets. Superior implementations exist in R already, however the `fastMCD` package represents my attempt to independently implement FAST-MCD and understand the underlying algorithm.


## Installation

```r
devtools::install_github("frankp-0/fastMCD", build_vignettes = T)
```

## Usage
The following examples demonstrate the usage of `fastMCD` for small (350 observations) and large (5000 observations) datasets.

```r
smallX <- MASS::mvrnorm(n = 350, mu = rep(0, 3), Sigma = diag(rep(1, 3)))
fastMCD(smallX)
```

```r
bigX <- MASS::mvrnorm(n = 5000, mu = rep(0, 3), Sigma = diag(rep(1, 3)))
fastMCD(bigX)
```

Users may optionally specify their own argument for h = # of observations to subsample.


```r
X <- MASS::mvrnorm(n = 350, mu = rep(0, 3), Sigma = diag(rep(1, 3)))
fastMCD(bigX, h = 500)
```
