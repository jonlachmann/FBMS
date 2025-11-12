[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/last-commit/jonlachmann/FBMS.svg)](https://github.com/jonlachmann/FBMS/commits/master)
[![](https://img.shields.io/github/languages/code-size/jonlachmann/FBMS.svg)](https://github.com/jonlachmann/FBMS)
[![R build status](https://github.com/jonlachmann/FBMS/workflows/R-CMD-check/badge.svg)](https://github.com/jonlachmann/FBMS/actions)
[![codecov](https://codecov.io/gh/jonlachmann/FBMS/branch/master/graph/badge.svg)](https://codecov.io/gh/jonlachmann/FBMS)
[![License: GPL](https://img.shields.io/badge/license-GPL-blue.svg)](https://cran.r-project.org/web/licenses/GPL)

# FBMS - Flexible Bayesian Model Selection

The `FBMS` package provides functions to estimate Bayesian Generalized nonlinear models (BGNLMs) through a Genetically Modified Mode Jumping MCMC algorithm.

# Installation and getting started
To install and load the development version of the package, just run
```
library(devtools)
install_github("jonlachmann/FBMS", force=T, build_vignettes=T)
library(FBMS)
```
With the package loaded, a vignette that shows how to run the package is available by running
```
vignette("FBMS-guide")
```
