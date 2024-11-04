# Introduction for R package TestSQR
We developed an hypothesis testing procedure for high-dimensional quantile regression model with high-dimensional nuisance parameters.
# Getting Started
These instructions will give you a toy example for implementing the package.
## Prerequisites
```
install.packages("devtools")
```
## Install TestSQR

```
library("devtools")
devtools::install_github("Bowen-stat/TestSQR")
```
## Toy example 
```
rm(list = ls())
set.seed(123)
library(TestSQR)