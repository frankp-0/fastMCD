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
The following examples demonstrate the usage of `fastMCD` for small (350 observations) and large (5000 observations) datasets which contain outliers. In both cases, the true means are 0 and the true covariance matrices are 3x3 identity matrices.

```r
library(fastMCD)
smallX <- MASS::mvrnorm(n = 350, mu = rep(0, 3), Sigma = diag(rep(1, 3)))
fastMCD(smallX)
## $center
## [1] -0.01491355  0.02881544  0.12785326
## 
## $cov
##             [,1]       [,2]        [,3]
## [1,]  0.47545610 -0.3179090 -0.07468779
## [2,] -0.31790899  1.7742387  0.19795384
## [3,] -0.07468779  0.1979538  1.49538074
```

```r
bigX <- MASS::mvrnorm(n = 350, mu = rep(0, 3), Sigma = diag(rep(1, 3)))
fastMCD(bigX)
## $center
## [1]  0.08968166 -0.07026145  0.18408677
## 
## $cov
##             [,1]        [,2]        [,3]
## [1,]  1.36256516  0.03485431 -0.03162824
## [2,]  0.03485431  0.73979040 -0.28927492
## [3,] -0.03162824 -0.28927492  1.17373954
```
