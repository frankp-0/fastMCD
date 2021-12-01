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
devtools::install_github("frankp-0/fastMCD")
```

## Usage
The following examples demonstrate the usage of `fastMCD` for small (350 observations) and large (5000 observations) datasets. In both cases, the true means are 0 and the true covariance matrices are 3x3 identity matrices.

```r
library(fastMCD)
data(smallX)
fastMCD(smallX)
## $center
## [1] -0.10337099 -0.09114819 -0.10676577
## 
## $cov
##               [,1]         [,2]        [,3]
## [1,]  0.8022050604 0.0008939482 -0.14231668
## [2,]  0.0008939482 0.8986171142  0.01015586
## [3,] -0.1423166796 0.0101558647  0.98102785
```

```r
data(bigX)
fastMCD(bigX)
## $center
## [1]  0.007787829 -0.021332101 -0.007643344
## 
## $cov
##             [,1]        [,2]         [,3]
## [1,]  0.97082365 0.028975299 -0.013889092
## [2,]  0.02897530 0.969469104  0.009680395
## [3,] -0.01388909 0.009680395  0.948956815
```
