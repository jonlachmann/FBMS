# Title     : Nonlinear functions for GMJMCMC
# Objective : Give an example library of nonlinear functions to use for GMJMCMC
# Created by: jonlachmann
# Created on: 2021-03-22

#' Sigmoid Function
#'
#' @param x The vector of values
#' @return The sigmoid of x
#'
#' @examples
#' sigmoid(2)
#' 
#'
#' @export sigmoid
sigmoid <- function(x)  1 / (1 + exp(-x))

#' Sine Function for Degrees
#'
#' @param x The vector of values in degrees
#' @return The sine of x
#'
#' @examples
#' sin_deg(0)
#'
#' @export sin_deg
sin_deg <- function(x) sin(x / 180 * pi)

#' Cosine Function for Degrees
#'
#' @param x The vector of values in degrees
#' @return The cosine of x
#'
#' @examples
#' cos_deg(0)
#'
#' @export cos_deg
cos_deg <- function(x) cos(x / 180 * pi)

#' Double Exponential Function
#'
#' @param x The vector of values
#' @return e^(-abs(x))
#'
#' @examples
#' exp_dbl(2)
#'
#' @export exp_dbl
exp_dbl <- function(x) exp(-abs(x))


#' Square Root Function
#'
#' @param x The vector of values
#' @return The square root of the absolute value of x
#'
#' @examples
#' sqroot(4)
#'
#' @export sqroot
sqroot <- function(x) abs(x)^(1/2)

#' Cube Root Function
#'
#' @param x The vector of values
#' @return The cube root of x
#'
#' @examples
#' troot(27)
#'
#' @export troot
troot <- function(x) abs(x)^(1/3)

#' Arcsinh Transform
#'
#' @param x The vector of values
#' @return arcsinh(x)
#'
#' @examples
#' arcsinh(2)
#'
#' @export arcsinh
arcsinh <- function(x) asinh(x)

#' p0 Polynomial Term
#'
#' @param x The vector of values
#' @return log(abs(x) + .Machine$double.eps)
#'
#' @examples
#' p0(2)
#'
#' @export p0
p0 <- function(x) log(abs(x)+.Machine$double.eps)

#' pm1 Polynomial Term
#'
#' @param x The vector of values
#' @return sign(x)*(abs(x)+.Machine$double.eps)^(-1)
#'
#' @examples
#' pm1(2)
#'
#' @export pm1
pm1 <- function(x) sign(x)*(abs(x)+.Machine$double.eps)^(-1)

#' pm2 Polynomial Term
#'
#' @param x The vector of values
#' @return sign(x)*(abs(x)+.Machine$double.eps)^(-2)
#'
#' @examples
#' pm2(2)
#'
#' @export pm2
pm2 <- function(x) sign(x)*(abs(x)+.Machine$double.eps)^(-2)

#' pm05 Polynomial Term
#'
#' @param x The vector of values
#' @return  (abs(x)+.Machine$double.eps)^(-0.5)
#'
#' @examples
#' pm05(2)
#'
#' @export pm05
pm05 <- function(x) (abs(x)+.Machine$double.eps)^(-0.5)

#' p05 Polynomial Term
#'
#' @param x The vector of values
#' @return (abs(x)+.Machine$double.eps)^(0.5)
#'
#' @examples
#' p05(2)
#'
#' @export p05
p05 <- function(x) (abs(x)+.Machine$double.eps)^(0.5)

#' p2 Polynomial Term
#'
#' @param x The vector of values
#' @return x^(2)
#'
#' @examples
#' p2(2)
#'
#' @export p2
p2 <- function(x) x^(2)

#' p3 Polynomial Term
#'
#' @param x The vector of values
#' @return x^(3)
#'
#' @examples
#' p3(2)
#'
#' @export p3
p3 <- function(x) x^(3)

#' p0p0 Polynomial Term
#'
#' @param x The vector of values
#' @return p0(x)*p0(x)
#'
#' @examples
#' p0p0(2)
#'
#' @export p0p0
p0p0 <- function(x) p0(x)*p0(x)

#' p0pm1 Polynomial Terms
#'
#' @param x The vector of values
#' @return p0(x)*(x+.Machine$double.eps)^(-1)
#'
#' @examples
#' p0pm1(2)
#'
#' @export p0pm1
p0pm1 <- function(x) p0(x)*(x+.Machine$double.eps)^(-1)

#' p0pm2 Polynomial Term
#'
#' @param x The vector of values
#' @return p0(x)*sign(x)*(abs(x)+.Machine$double.eps)^(-2)
#'
#' @examples
#' p0pm2(2)
#'
#' @export p0pm2
p0pm2 <- function(x) p0(x)*sign(x)*(abs(x)+.Machine$double.eps)^(-2)

#' p0pm05 Polynomial Term
#'
#' @param x The vector of values
#' @return p0(x)*sign(x)*(abs(x)+.Machine$double.eps)^(-0.5)
#'
#' @examples
#' p0pm05(2)
#'
#' @export p0pm05
p0pm05 <- function(x) p0(x)*(abs(x)+.Machine$double.eps)^(-0.5)

#' p0p05 Polynomial Term
#'
#' @param x The vector of values
#' @return p0(x)*(abs(x)+.Machine$double.eps)^(0.5)
#'
#' @examples
#' p0p05(2)
#'
#' @export p0p05
p0p05 <- function(x) p0(x)*(abs(x)+.Machine$double.eps)^(0.5)

#' p0p1 Polynomial Term
#'
#' @param x The vector of values
#' @return  p0(x)*x
#'
#' @examples
#' p0p1(2)
#'
#' @export p0p1
p0p1 <- function(x) p0(x)*x

#' p0p2 Polynomial Term
#'
#' @param x The vector of values
#' @return p0(x)*x^(2)
#'
#' @examples
#' p0p2(2)
#'
#' @export p0p2
p0p2 <- function(x) p0(x)*x^(2)

#' p0p3 Polynomial Term
#'
#' @param x The vector of values
#' @return p0(x)*x^(3)
#'
#' @examples
#' p0p3(2)
#'
#' @export p0p3
p0p3 <- function(x) p0(x)*x^(3)


#' ReLU Function
#'
#' @param x The vector of values
#' @return max(x,0)
#'
#' @examples
#' relu(2)
#'
#' @export relu
relu <- function(x) max(x,0)


#' Negative ReLU Function
#'
#' @param x The vector of values
#' @return max(-x,0)
#'
#' @examples
#' nrelu(2)
#'
#' @export nrelu
nrelu <- function(x) max(-x,0)

#' GELU Function
#'
#' @param x The vector of values
#' @return x*pnorm(x)
#'
#' @examples
#' gelu(2)
#'
#' @export gelu
gelu <- function(x)x *pnorm(x)


#' Negative GELU Function
#'
#' @param x The vector of values
#' @return -x*pnorm(-x)
#'
#' @examples
#' ngelu(2)
#'
#' @export ngelu
ngelu <- function(x) -x*pnorm(-x)

#' Erf Function
#'
#' @param x The vector of values
#' @return 2 * pnorm(x * sqrt(2)) - 1
#'
#' @examples
#' erf(2)
#'
#' @export erf
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


#' Heavy Side Function
#'
#' @param x The vector of values
#' @return as.integer(x>0)
#'
#' @examples
#' hs(2)
#'
#' @export hs
hs <- function(x) as.integer(x>0)


#' Negative Heavy Side Function
#'
#' @param x The vector of values
#' @return as.integer(x<0)
#'
#' @examples
#' nhs(2)
#'
#' @export nhs
nhs <- function(x) as.integer(x<0)

#' Not x
#'
#' @param x The vector of binary values
#' @return 1-x
#'
#' @examples
#' not(TRUE)
#'
#' @export not
not <- function(x) (1-x)
