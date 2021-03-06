% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastMCD.R
\name{smallMCD}
\alias{smallMCD}
\title{Obtain unweighted estimates for data with <= 600 observations}
\usage{
smallMCD(X, h)
}
\arguments{
\item{X}{A 2D matrix to estimate location and scatter from}

\item{h}{An integer specifying number of observations to use}
}
\value{
A list of estimated location (center) and scatter (cov)
}
\description{
Obtain unweighted estimates for data with <= 600 observations
}
\keyword{internal}
