# Title     : Nonlinear functions for GMJMCMC
# Objective : Give an example library of nonlinear functions to use for GMJMCMC
# Created by: jonlachmann
# Created on: 2021-03-22

#' Sigmoid function
#'
#' @param x The vector of values
#' @return The sigmoid of x
#'
#' @export sigmoid
sigmoid <- function(x) exp(-x)

#' Sine function for degrees
#'
#' @param x The vector of values in degrees
#' @return The sine of x
#'
#' @export sini
sini <- function(x) sin(x/180*pi)

#' Cosine function for degrees
#'
#' @param x The vector of values in degrees
#' @return The cosine of x
#'
#' @export cosi
cosi <- function(x) cos(x/180*pi)

#' Exponential function of absolute values
#'
#' @param x The vector of values
#' @return e^(-abs(x))
#'
#' @export expi
expi <- function(x) exp(-abs(x))

#' Log function of absolute values
#'
#' @param x The vector of values
#' @return log(abs(x))
#'
#' @export logi
logi <- function(x) log(abs(x)+.Machine$double.eps)

#' Square root function
#'
#' @param x The vector of values
#' @return The square root of the absolute value of x
#'
#' @export sqroot
sqroot <- function(x) abs(x)^(1/2)

#' Cube root function
#'
#' @param x The vector of values
#' @return The cube root of x
#'
#' @export troot
troot <- function(x) abs(x)^(1/3)

#' To the 2.3  power function
#'
#' @param x The vector of values
#' @return x^2.3
#'
#' @export troot
to23 <- function(x) abs(x)^(2.3)

#' To the 7/2  power function
#'
#' @param x The vector of values
#' @return x^(7/2)
#'
#' @export troot
to72 <- function(x) abs(x)^(7/2)

#' Gaussian function
#'
#' @param x The vector of values
#' @return e^(-x^2)
#'
#' @export gauss
gauss <- function(x) exp(-x*x)

#' To 2.5 power
#'
#' @param x The vector of values
#' @return x^(2.5)
#'
#' @export to25
to25 <- function(x)abs(x)^(2.5)

#' To 3 power
#'
#' @param x The vector of values
#' @return x^(3)
#'
#' @export to35
to3 <- function(x)abs(x)^(3)

#' To 3.5 power
#'
#' @param x The vector of values
#' @return x^(3.5)
#'
#' @export to35
to35 <- function(x)abs(x)^(3.5)