% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{soft}
\alias{soft}
\title{Element-wise Soft Thresholding Function for Matrix}
\usage{
soft(X, lambda = 0)
}
\arguments{
\item{\code{X}, }{a Matrix}

\item{\code{lambda}, }{a scalar}
}
\value{
Matrix after threholding
}
\description{
Soft-thresholding function for a matrix which is the solution of
\deqn{\argmin_Y ||Y-X||_2^2/2+\lambda ||Y||_1.}
}
