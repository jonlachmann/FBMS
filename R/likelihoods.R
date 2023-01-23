# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24

#' Log likelihood function for logistic regression with a prior p(m)=sum(total_width)
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @export logistic.loglik
logistic.loglik <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=binomial())})
  ret <- (-(mod$deviance -2*log(params$r)*sum(complex$oc)))/2
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log likelihood function for logistic regression for alpha calculation
#' This function is just the bare likelihood function
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a[i].
#'
#' @export logistic.loglik.alpha
logistic.loglik.alpha <- function (a, data, mu_func) {
  m <- 1/(1+exp(-eval(parse(text=mu_func))))
  -sum((data[,1] * log(m) + (1-data[,1]) * log(1 - m)))
}

#' Log likelihood function for gaussian regression with a prior p(m)=r*sum(total_width).
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @export gaussian.loglik
gaussian.loglik <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=gaussian())})
  ret <- (-(mod$deviance + 2*log(length(y))*(mod$rank-1) - 2*log(params$r)*(sum(complex$oc))))/2
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log likelihood function for gaussian regression with a prior p(m)=r*sum(total_width), using subsampling.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @export gaussian.loglik.bic.irlssgd
gaussian.loglik.bic.irlssgd <- function (y, x, model, complex, params) {
  mod <- irls.sgd(as.matrix(x[,model]), y, gaussian(),
            irls.control=list(subs=params$subs, maxit=20, tol=1e-7, cooling = c(1,0.9,0.75), expl = c(3,1.5,1)),
            sgd.control=list(subs=params$subs, maxit=250, alpha=0.001, decay=0.99, histfreq=10))
  ret <- (-(mod$deviance + 2*log(length(y))*(mod$rank-1) - 2*log(params$r)*(sum(complex$oc))))/2
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log likelihood function for gaussian regression for alpha calculation
#' This function is just the bare likelihood function
#' Note that it only gives a proportional value and is equivalent to least squares
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a[i].
#'
#' @export gaussian.loglik.alpha
gaussian.loglik.alpha <- function (a, data, mu_func) {
  m <- eval(parse(text=mu_func))
  sum((data[,1]-m)^2)
}

#' Log likelihood function for linear regression using Zellners g-prior
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @export linear.g.prior.loglik
linear.g.prior.loglik <- function (y, x, model, complex, params) {
  out <- lm.fit(as.matrix(x[,model]),y)
  rsquared <- 1-sum(var(out$residuals))/sum(var(y))
  p <- out$rank
  n <- nrow(x)
  logmarglik <- 0.5*(log(1+params$g)*(n-p) - log(1+params$g*(1-rsquared))*(n-1))*(p!=1)
  return(logmarglik)
}

#' Log likelihood function for a VAR model
#'
#' @param y A matrix containing the dependent variables
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @export logistic.loglik
var.loglik <- function (y, x, model, complex, params) {
  # Rearrange the data as it is multivariate
  y <- cbind(y, x[, seq_len(params$K-1)])
  x <- x[, -seq_len(params$K-1)]
  r <- 20/223
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=binomial())})
  ret <- (-(mod$deviance -2*log(r)*sum(complex$oc)))/2
  return(ret)
}