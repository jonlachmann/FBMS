# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24

#' Log likelihood function for glm regression with parameter priors from BAS package
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user, important to specify the tuning parameters of beta priors and family that BAS uses in glm models
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' glm.logpost.bas(as.integer(rnorm(100) > 0),cbind(1,matrix(rnorm(100))),c(TRUE,TRUE),list(oc = 1))
#'
#' @importFrom BAS uniform Jeffreys g.prior
#' @importFrom stats poisson Gamma glm.control
#' @export glm.logpost.bas
glm.logpost.bas2 <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5), family = "binomial", prior_beta = Jeffreys(), laplace = FALSE)) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r =  1/dim(x)[1], family = "binomial", prior_beta = g.prior(max(dim(x)[1], sum(model)-1)), laplace = FALSE)
  p <- sum(model) - 1
  if(p==0)
  {
    probinit <- as.numeric(c(1,0.99))
    model[2] <- T
  }else{
    probinit <- as.numeric(c(1,rep(0.99,p)))
  }

  mod<-NULL

  
  tryCatch({
  if(mlpost_params$family == "binomial")
    suppressWarnings({
      mod <- .Call(BAS:::C_glm_deterministic,
                   y = as.numeric(y), X = as.matrix(x[,model]),
                                  Roffset = as.numeric(rep(0, length(y))),
                                  Rweights = as.numeric(rep(1, length(y))),
                                  Rprobinit = probinit,
                                  Rmodeldim = as.integer(rep(0,ifelse(p==0,2,1))),
                                                         modelprior = uniform(),
                                                         betaprior = mlpost_params$prior_beta,
                                                         family = binomial(),
                                                         Rcontrol = glm.control(),
                                                         Rlaplace =  as.integer(mlpost_params$laplace))
    })
  else if(mlpost_params$family == "poisson")
    suppressWarnings({
      mod <- .Call(BAS:::C_glm_deterministic,
                   y = as.numeric(y),  X = as.matrix(x[,model]),
                                  Roffset = as.numeric(rep(0, length(y))),
                                  Rweights = as.numeric(rep(1, length(y))),
                                  Rprobinit = probinit,
                                  Rmodeldim = as.integer(rep(0,ifelse(p==0,2,1))),
                                                         modelprior = uniform(),
                                                         betaprior = mlpost_params$prior_beta,
                                                         family = poisson(),
                                                         Rcontrol = glm.control(),
                                                         Rlaplace =  as.integer(mlpost_params$laplace))
    })
  else{
    suppressWarnings({
      mod <- .Call(BAS:::C_glm_deterministic,
                   y = as.numeric(y), X = as.matrix(x[,model]),
                                  Roffset = as.numeric(rep(0, length(y))),
                                  Rweights = as.numeric(rep(1, length(y))),
                                  Rprobinit = probinit,
                                  Rmodeldim = as.integer(rep(0,ifelse(p==0,2,1))),
                                                         modelprior = uniform(),
                                                         betaprior = mlpost_params$prior_beta,
                                                         family = Gamma(),
                                                         Rcontrol = glm.control(),
                                                         Rlaplace =  as.integer(mlpost_params$laplace))
    })}
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })

  if(length(mod)==0) {
    return(list(crit = -.Machine$double.xmax + log(mlpost_params$r * sum(complex$oc)), coefs = rep(0, p+1)))
  }

  if(p == 0)
  {
    ret <- mod$logmarg[2] + log(mlpost_params$r) * sum(complex$oc)
    return(list(crit=ret, coefs=mod$mle[[2]]))
  }
  ret <- mod$logmarg + log(mlpost_params$r) * sum(complex$oc)
  return(list(crit=ret, coefs=mod$mle[[1]]))
}


#' Log likelihood function for Gaussian regression with parameter priors from BAS package
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user, important to specify the tuning parameters of beta priors where the corresponding integers as prior_beta must be provided "g-prior" = 0, "hyper-g" = 1, "EB-local" = 2, "BIC" = 3, "ZS-null" = 4, "ZS-full" = 5, "hyper-g-laplace" = 6, "AIC" = 7, "EB-global" = 2, "hyper-g-n" = 8, "JZS" = 9 and in Gaussian models
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' lm.logpost.bas(rnorm(100), cbind(1,matrix(rnorm(100))), c(TRUE,TRUE), list(oc = 1))
#'
#'
#' @export lm.logpost.bas
lm.logpost.bas2 <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5), prior_beta = "g-prior", alpha = 4)) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r =  1/dim(x)[1], prior_beta = 0, alpha = max(dim(x)[1], sum(model)^2))


  p <- sum(model) - 1
  if(p==0)
  {
    probinit <- as.numeric(c(1,0.99))
    model[2] <- T
  }else{
    probinit <- as.numeric(c(1,rep(0.99,p)))
  }

  mod<-NULL

  tryCatch({
      suppressWarnings({
        mod <- .Call(BAS:::C_deterministic,
                     y = y, X = as.matrix(x[,model]),
                     as.numeric(rep(1, length(y))),
                     probinit,
                     as.integer(rep(0,ifelse(p==0,2,1))),
                     incint = as.integer(F),
                     alpha = ifelse(length(mlpost_params$alpha)>0, as.numeric(mlpost_params$alpha), NULL),
                     method = as.integer(mlpost_params$prior_beta),
                     modelprior = uniform(),
                     Rpivot = TRUE,
                     Rtol = 1e-7)
      })
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })

  if(length(mod)==0) {
    return(list(crit = -.Machine$double.xmax + log(mlpost_params$r * sum(complex$oc)), coefs = rep(0, p+1)))
  }

  if(p == 0)
  {
    ret <- mod$logmarg[2] + log(mlpost_params$r) * sum(complex$oc)
    return(list(crit=ret, coefs=mod$mle[[2]]))
  }
  ret <- mod$logmarg + log(mlpost_params$r) * sum(complex$oc)
  return(list(crit=ret, coefs=mod$mle[[1]]))
}


#' Log likelihood function for logistic regression with a Jeffreys parameter prior and BIC approximations of the posterior
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' logistic.loglik2(as.integer(rnorm(100) > 0), matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export logistic.loglik2
logistic.loglik2 <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5))) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r = 1/dim(x)[1])
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = binomial())})
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(mlpost_params$r) * sum(complex$oc))) / 2
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log likelihood function for glm regression with a Jeffreys parameter prior and BIC approximations of the posterior
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user, family must be specified
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' glm.loglik(abs(rnorm(100))+1, matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export glm.loglik
glm.loglik2 <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5), family = "Gamma")) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r = 1/dim(x)[1])

  if(mlpost_params$family == "binomial")
  {
    fam = binomial()
  }else if(mlpost_params$family == "poisson"){
    fam = poisson()
  }else
  {
    fam = Gamma()
  }

  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = fam)})
  
  if (length(mod) == 0 || is.nan(mod$deviance)) {
    return(list(crit = -.Machine$double.xmax + log_prior(mlpost_params, complex), coefs = rep(0, sum(model))))
  }
  
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(mlpost_params$r) * sum(complex$oc))) / 2
  return(list(crit=ret, coefs=mod$coefficients))
}


#' Log likelihood function for gaussian regression with a Jeffreys prior and BIC approximation of MLIK with both known and unknown variance of the responses
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' gaussian.loglik(rnorm(100), matrix(rnorm(100)), TRUE, list(oc = 1), NULL)
#'
#'
#' @export gaussian.loglik
gaussian.loglik2 <- function (y, x, model, complex, mlpost_params) {
  if(length(mlpost_params)==0)
    mlpost_params <- list()
  if (length(mlpost_params$r) == 0)
    mlpost_params$r <- 1/dim(x)[1]
  if(length(mlpost_params$var) == 0)
    mlpost_params$var <- 1
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())})

  if(mlpost_params$var == "unknown")
    ret <- (-(mod$aic + (log(length(y))-2) * (mod$rank) - 2 * log(mlpost_params$r) * (sum(complex$oc)))) / 2
  else
    ret <- (-(mod$deviance/mlpost_params$var + log(length(y)) * (mod$rank - 1) - 2 * log_prior(mlpost_params, complex))) / 2

  return(list(crit=ret, coefs=mod$coefficients))
}


#' Log likelihood function for linear regression using Zellners g-prior
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' gaussian.loglik2.g(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @export gaussian.loglik2.g
gaussian.loglik2.g <- function (y, x, model, complex, mlpost_params = NULL)
{
  if(length(mlpost_params)==0)
    mlpost_params <- list()
  if (length(mlpost_params$r) == 0)
    mlpost_params$r <- 1/dim(x)[1]
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  # Calculate R-squared
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  Rsquare <- 1 - (RSS / TSS)

  if (length(mlpost_params$r) == 0 || length(mlpost_params$g) == 0)
  {
    mlpost_params$r <- 1/dim(x)[1]
    mlpost_params$g <- max(mod$rank^2, length(y))
  }

  # logarithm of marginal likelihood
  mloglik <- 0.5*(log(1.0 + mlpost_params$g) * (dim(x)[1] - mod$rank)  - log(1.0 + mlpost_params$g * (1.0 - Rsquare)) * (dim(x)[1]  - 1))*(mod$rank!=1)

  # logarithm of model prior
   # default value or parameter r
  lp <- log_prior(mlpost_params, complex)

  return(list(crit = mloglik + lp, coefs = mod$coefficients))
}


#' Log likelihood function for Gaussian regression with parameter priors from BAS package
#'
#' This function computes the marginal likelihood of a Gaussian regression model under different priors.
#'
#' @param y A numeric vector containing the dependent variable.
#' @param x A matrix containing the independent variables, including an intercept column.
#' @param model A logical vector indicating which variables to include in the model.
#' @param complex A list containing complexity measures for the features.
#' @param mlpost_params A list of parameters for the log likelihood, specifying the tuning parameters of beta priors.
#'
#' @return A list with elements:
#'   \item{crit}{Log marginal likelihood combined with the log prior.}
#'   \item{coefs}{Posterior mode of the coefficients.}
#'
#' @examples
#' gaussian_tcch_log_likelihood2(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @importFrom BAS phi1 hypergeometric1F1 hypergeometric2F1
#' @importFrom tolerance F1
#' @export
gaussian_tcch_log_likelihood2 <- function(y, x, model, complex, mlpost_params = list(r = exp(-0.5), prior_beta = "intrinsic")) {

  # Fit the linear model using fastglm
  fitted_model <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  log_likelihood <- -(fitted_model$aic  -2 * (fitted_model$rank))/2
  # Compute R-squared manually
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(fitted_model$residuals^2)
  R2_M <- 1 - (RSS / TSS)

  p_M <- fitted_model$rank
  n <- length(y)

  # Switch-like structure to assign hyperparameters based on prior
  if (mlpost_params$prior_beta[[1]] == "CH") {
    # CH prior: b and s should be user-specified, with defaults if not provided
    a <- ifelse(!is.null(mlpost_params$prior_beta$a), mlpost_params$prior_beta$a, 1)  # Default to 1 if not specified
    b <- ifelse(!is.null(mlpost_params$prior_beta$b), mlpost_params$prior_beta$b, 2)  # Default to 1 if not specified
    r <- 0
    s <- ifelse(!is.null(mlpost_params$prior_beta$s), mlpost_params$prior_beta$s, 1)  # Default to 1 if not specified
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "hyper-g") {
    a <- 1
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "uniform") {
    a <- 2
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "Jeffreys") {
    a <- 0.0001
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$prior_beta[[1]] == "beta.prime") {
    a <- 1/2
    b <- n - p_M - 1.5
    r <- 0
    s <- 0
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "benchmark") {
    a <- 0.02
    b <- 0.02 * max(n, p_M^2)
    r <- 0
    s <- 0
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "TG") {

    a <- 2 * ifelse(!is.null(mlpost_params$prior_beta$a), mlpost_params$prior_beta$a, 1)
    b <- 2
    r <- 0
    s <- 2 * ifelse(!is.null(mlpost_params$prior_beta$s), mlpost_params$prior_beta$s, 1)
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "ZS-adapted") {
    a <- 1
    b <- 2
    r <- 0
    s <- n + 3
    v <- 1
    k <- 1
  } else if (mlpost_params$prior_beta[[1]] == "robust") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- (n + 1) / (p_M + 1)
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "hyper-g-n") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- 1
    k <- 1

  } else if (mlpost_params$prior_beta[[1]] == "intrinsic") {
    a <- 1
    b <- 1
    r <- 1
    s <- 0
    v <- (n + p_M + 1) / (p_M + 1)
    k <- (n + p_M + 1) / n

  }else if (mlpost_params$prior_beta[[1]] == "tCCH") {
    a <- mlpost_params$prior_beta$a
    b <- mlpost_params$prior_beta$b
    r <- mlpost_params$prior_beta$rho
    s <- mlpost_params$prior_beta$s
    v <- mlpost_params$prior_beta$v
    k <- mlpost_params$prior_beta$k
  }else {
    stop("Unknown prior name: ", mlpost_params$prior_beta)
  }

  #
  if (!is.null(r) & r == 0) {

    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- phi1(b / 2, (n - 1) / 2, (a + b + p_M) / 2, s / (2 * v), min(0.8,R2_M/(v - (v - 1) * R2_M),log = T))

    if(R2_M/(v - (v - 1) * R2_M)>0.8)
    {
      warning("Infinite marginal log likelihood! phi1 last argument reduced to 0.8. Use a different prior_beta (Robust, Hyper-g/n, Intrinsic, or g-prior)")
    }

    term3 <- lbeta(a / 2, b / 2)
    term4 <- hypergeometric1F1(b / 2, (a + b) / 2, s / (2 * v),log = T)
    marginal_likelihood <- log_likelihood + (term1) + (term2) - (p_M / 2) * log(v) - ((n - 1) / 2)*log(1 - (1 - 1 / v) * R2_M) - (term3) - (term4)
  } else if (!is.null(s) & s == 0) {
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- hypergeometric2F1(r, b / 2, (a + b) / 2, 1 - k,log = T)
    term3 <- F1((a + p_M) / 2, (a + b + p_M + 1 - n - 2 * r) / 2, (n - 1) / 2, (a + b + p_M) / 2, 1 - k, 1 - k - (R2_M^2 * k) / ((1 - R2_M) * v))
    marginal_likelihood <- log_likelihood + (a+p_M-2*r)/2*log(k) + (term1) - (term2) - (term3) - (p_M / 2) * log(v) - log(1 - R2_M) * ((n - 1) / 2) - lbeta(a / 2, b / 2)

  } else {
    stop("Invalid inputs: either r = 0 or s = 0 must be specified.")
  }

  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r

  lp <- log_prior(mlpost_params, complex)

  return(list(crit = marginal_likelihood + lp, coefs = fitted_model$coefficients))
}



#' Log likelihood function for logistic regression with an approximate Laplace approximations used
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' logistic.loglik2.ala(as.integer(rnorm(100) > 0), matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export logistic.loglik2.ala
logistic.loglik2.ala <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5))) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r = 1/dim(x)[1])
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = binomial(),maxit = 1)})
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) -2 * log(mlpost_params$r) * sum(complex$oc))) / 2
  return(list(crit=ret, coefs=mod$coefficients))
}



#' Log likelihood function for logistic regression for alpha calculation
#' This function is just the bare likelihood function
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a\[i\].
#'
#' @return A numeric with the log likelihood.
#'
#' @export logistic.loglik2.alpha
logistic.loglik2.alpha <- function (a, data, mu_func) {
  m <- 1 / (1 + exp(-eval(parse(text = mu_func))))
  -sum((data[,1] * log(m) + (1 - data[, 1]) * log(1 - m)))
}


#' Log likelihood function for gaussian regression for alpha calculation
#' This function is just the bare likelihood function
#' Note that it only gives a proportional value and is equivalent to least squares
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a\[i\].
#'
#' @return A numeric with the log likelihood.
#' @examples
#'\dontrun{
#'gaussian.loglik2.alpha(my_alpha,my_data,my_mu)
#'}
#' @export gaussian.loglik2.alpha
gaussian.loglik2.alpha <- function (a, data, mu_func) {
  m <- eval(parse(text=mu_func))
  sum((data[,1]-m)^2)
}


#' Log model prior function
#' @param mlpost_params list of passed parameters of the likelihood in GMJMCMC
#' @param complex list of complexity measures of the features included into the model
#'
#' @return A numeric with the log model prior.
#'
#' @examples
#' log_prior(mlpost_params = list(r=2), complex = list(oc = 2))
#'
#' @export log_prior
log_prior <- function (mlpost_params, complex) {
  pl <- log(mlpost_params$r) * (sum(complex$oc))
  return(pl)
}


#' Master Log Marginal Likelihood Function
#'
#' This function serves as a unified interface to compute the log marginal likelihood
#' for different regression models and priors by calling specific log likelihood functions.
#'
#' @param y A numeric vector containing the dependent variable.
#' @param x A matrix containing the precalculated features (independent variables).
#' @param model A logical vector indicating which variables to include in the model.
#' @param complex A list of complexity measures for the features.
#' @param mlpost_params A list of parameters controlling the model family, prior, and tuning parameters.
#'   Key elements include:
#'   - family: "binomial", "poisson", "gamma" (all three referred to as GLM below), or "gaussian" (default: "gaussian")
#'   - prior_beta: Type of prior as a string (default: "g-prior"). Possible values include:
#'     - "beta.prime": Beta-prime prior (GLM/Gaussian, no additional args)
#'     - "CH": Compound Hypergeometric prior (GLM/Gaussian, requires `a`, `b`, optionally `s`)
#'     - "EB-local": Empirical Bayes local prior (GLM/Gaussian, requires `a` for Gaussian)
#'     - "EB-global": Empirical Bayes local prior (Gaussian, requires `a` for Gaussian)
#'     - "g-prior": Zellner's g-prior (GLM/Gaussian, requires `g`)
#'     - "hyper-g": Hyper-g prior (GLM/Gaussian, requires `a`)
#'     - "hyper-g-n": Hyper-g/n prior (GLM/Gaussian, requires `a`)
#'     - "tCCH": Truncated Compound Hypergeometric prior (GLM/Gaussian, requires `a`, `b`, `s`, `rho`, `v`, `k`)
#'     - "intrinsic": Intrinsic prior (GLM/Gaussian, no additional args)
#'     - "TG": Truncated Gamma prior (GLM/Gamma, requires `a`, `s`)
#'     - "Jeffreys": Jeffreys prior (GLM/Gaussian, no additional args)
#'     - "uniform": Uniform prior (GLM/Gaussian, no additional args)
#'     - "benchmark": Benchmark prior (Gaussian/GLM, no additional args)
#'     - "ZS-adapted": Zellner-Siow adapted prior (Gaussian TCCH, no additional args)
#'     - "robust": Robust prior (Gaussian/GLM, no additional args)
#'     - "Jeffreys-BIC": Jeffreys prior with BIC approximation of marginal likelihood (Gaussian/GLM)
#'     - "ZS-null": Zellner-Siow null prior (Gaussian, requires `a`)
#'     - "ZS-full": Zellner-Siow full prior (Gaussian, requires `a`)
#'     - "hyper-g-laplace": Hyper-g Laplace prior (Gaussian, requires `a`)
#'     - "AIC": AIC prior from BAS (Gaussian, requires penalty `a`)
#'     - "BIC": BIC prior from BAS (Gaussian/GLM)
#'     - "JZS": Jeffreys-Zellner-Siow prior (Gaussian, requires `a`)
#'   - r: Model complexity penalty (default: 1/n)
#'   - g: Tuning parameter for g-prior (default: max(n, p^2))
#'   - a, b, s, v, rho, k: Hyperparameters for various priors
#'   - n: Sample size for some priors (default: length(y))
#'   - var: Variance assumption for Gaussian models ("known" or "unknown", default: "unknown")
#'   - laplace: Logical for Laplace approximation in GLM only (default: FALSE)
#'
#' @return A list with elements:
#'   \item{crit}{Log marginal likelihood combined with the log prior.}
#'   \item{coefs}{Posterior mode of the coefficients.}
#'
#' @examples
#' fbms.mlik.master(rnorm(100), matrix(rnorm(100)), c(TRUE,TRUE), list(oc = 1))
#'
#' @importFrom BAS robust beta.prime bic.prior CCH EB.local g.prior hyper.g hyper.g.n tCCH intrinsic TG Jeffreys uniform
#' @export
fbms.mlik.master_old <- function(y, x, model, complex, mlpost_params = list(family = "gaussian", prior_beta = "g-prior", r = exp(-0.5))) {
  # Extract dimensions
  n <- length(y)
  p <- sum(model) - 1  # Number of predictors excluding intercept

  # Set default parameters if not fully specified
  if (is.null(mlpost_params$family)) mlpost_params$family <- "gaussian"
  if (is.null(mlpost_params$prior_beta)) mlpost_params$prior_beta <- "g-prior"
  if (is.null(mlpost_params$g)) mlpost_params$g <- max(p^2, n)
  if (is.null(mlpost_params$n)) mlpost_params$n <- n
  if (is.null(mlpost_params$r)) mlpost_params$r <- 1/n

  # Ensure complex has oc if not provided, ignore by default
  if (is.null(complex$oc)) complex$oc <- 0

  # Homogenize and prepare mlpost_params for nested calls
  params_master <- mlpost_params
  params_nested <- list(r = params_master$r)

  # Define valid priors for each family
  #glm_only_priors <- c("CCH", "tCCH", "TG")
  glm_and_gaussian_priors <- c("CH", "tCCH", "TG","beta.prime", "EB-local", "g-prior", "hyper-g", "hyper-g-n",
                               "intrinsic", "ZS-adapted", "Jeffreys", "uniform", "benchmark", "robust", "Jeffreys-BIC")
  gaussian_only_priors <- c("ZS-null", "ZS-full", "hyper-g-laplace","BIC", "AIC", "JZS","EB-global")

  #review a bit
  gaussian_not_robust <-  c("CH", "tCCH", "ZS-adapted", "TG","beta.prime", "benchmark","Jeffreys")
  gaussian_robust <- c("g-prior", "hyper-g", "EB-local","BIC", "Jeffreys-BIC", "ZS-null", "ZS-full", "hyper-g-laplace",
                       "AIC",  "hyper-g-n",  "JZS")
  gaussian_tcch <- c("CH", "tCCH", "TG","beta.prime", "intrinsic", "ZS-adapted", "uniform","Jeffreys", "benchmark", "robust")
  gaussian_bas <- c("g-prior", "hyper-g", "EB-local","ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS")

  all_priors <- c(glm_and_gaussian_priors, gaussian_only_priors)
  #browser()
  # Validate prior_beta
  if (!params_master$prior_beta %in% all_priors) {
    stop(sprintf("Prior '%s' is not supported. Supported priors: %s",
                 params_master$prior_beta, paste(all_priors, collapse = ", ")))
  }

  # Decision logic based on family and prior_beta
  if (params_master$family %in% c("binomial", "poisson", "gamma")) {
    if (params_master$prior_beta %in% gaussian_only_priors) {
      stop(sprintf("Prior '%s' is not supported for GLM family '%s'. Supported GLM priors: %s",
                   params_master$prior_beta, params_master$family,
                   paste(c(glm_and_gaussian_priors), collapse = ", ")))
    }

    params_nested$family <- params_master$family
    if (is.null(params_master$laplace)) params_nested$laplace <- FALSE else params_nested$laplace <- params_master$laplace

    #if(params_nested$laplace)
    #  print("Laplace approximations will be used")

    if (params_master$prior_beta == "Jeffreys-BIC") {
      if(params_nested$family == "binomial")
        result <- logistic.loglik2(y, x, model, complex, params_nested)
      else if(params_nested$family%in% c("poisson", "gamma"))
        result <- glm.loglik2(y, x, model, complex, params_nested)

    } else {
      params_nested$prior_beta <- switch(
        params_master$prior_beta,
        "beta.prime" = beta.prime(n = n),
        "CH" = CCH(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                    beta = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                    s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s),
        "EB-local" = EB.local(),
        "g-prior" = g.prior(g = params_master$g),
        "hyper-g" = hyper.g(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a),
        "tCCH" = tCCH(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                      beta = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                      s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s,
                      r = if (is.null(params_master$rho)) stop("rho must be provided") else params_master$rho,
                      v = if (is.null(params_master$v)) stop("v must be provided") else params_master$v,
                      theta = if (is.null(params_master$k)) stop("k must be provided") else params_master$k),
        "intrinsic" = intrinsic(n = params_master$n),
        "TG" = TG(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a),
        "Jeffreys" = Jeffreys(),
        "uniform" = tCCH(alpha = 2,
                         beta = 2,
                         s = 0,
                         r = 0,
                         v = 1,
                         theta = 1),
        "benchmark" =  tCCH(alpha = 0.02,
                            beta = 0.02*max(n,p^2),
                            s = 0,
                            r = 0,
                            v = 1,
                            theta = 1),
        "ZS-adapted" =  tCCH(alpha = 1,
                            beta = 2,
                            s = n + 3,
                            r = 0,
                            v = 1,
                            theta = 1),
        "TG" = TG(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a),
        "robust" = robust(n = if (is.null(params_master$n)) as.numeric(n) else as.numeric(params_master$n)),
        "hyper-g-n" = hyper.g.n(alpha = if (is.null(params_master$a)) 3 else params_master$a,
                                n = params_master$n),
        "BIC" = bic.prior(n = if (is.null(params_master$n)) n else params_master$n),
        stop("Unrecognized prior_beta for GLM: ", params_master$prior_beta)
      )
      result <- glm.logpost.bas2(y, x, model, complex, params_nested)
    }
  } else if (params_master$family == "gaussian") {

    if (params_master$prior_beta %in% gaussian_not_robust) {
      warning(sprintf("Prior '%s' is not reccomended supported for Gaussian family '%s'. Can be unstable for strong signals (R^2 > 0.9). Reccomended priors under Gaussian family: %s",
                   params_master$prior_beta, params_master$family,
                   paste(gaussian_robust, collapse = ", ")))
    }

    params_nested$r <- params_master$r

    if (params_master$prior_beta %in% gaussian_tcch) {

      params_nested$prior_beta <- switch(
        params_master$prior_beta,
        "beta.prime" = list("beta.prime"),
        "CH" = list("CH",a = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                    b = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                    s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s),
        "tCCH" = list("tCCH", a = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                      b = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                      s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s,
                      rho = if (is.null(params_master$rho)) stop("rho must be provided") else params_master$rho,
                      v = if (is.null(params_master$v)) stop("v must be provided") else params_master$v,
                      k = if (is.null(params_master$k)) stop("k must be provided") else params_master$k),
        "intrinsic" = list("intrinsic"),
        "TG" = list("TG",a = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                    s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s),
        "Jeffreys" = list("Jeffreys"),
        "ZS-adapted" = list("ZS-adapted"),
        "benchmark" = list("benchmark"),
        "robust" = list("robust"),
        "uniform" = list("uniform"),
        stop("Unrecognized prior_beta for Gaussian GLM: ", params_master$prior_beta)
      )
      result <- gaussian_tcch_log_likelihood2(y, x, model, complex, params_nested)

    }else if (params_master$prior_beta == "Jeffreys-BIC") {
      if (is.null(params_master$var)) params_nested$var <- "unknown" else params_nested$var <- params_master$var
      result <- gaussian.loglik2(y, x, model, complex, params_nested)
    } else if (params_master$prior_beta %in% gaussian_bas) {

      params_nested$prior_beta  <- switch(
        params_master$prior_beta,
        "g-prior" = 0,
        "hyper-g" = 1,
        "EB-local" = 2,
        "BIC" = 3,
        "ZS-null" = 4,
        "ZS-full" = 5,
        "hyper-g-laplace" = 6,
        "AIC" = 7,
        "EB-global" = 2,
        "hyper-g-n" = 8,
        "JZS" = 9
      )
      if(params_master$prior_beta == "g-prior")
      {
        if (!is.null(params_master$g)) params_nested$g <- params_master$g else stop("g must be provided")
        result <- gaussian.loglik2.g(y, x, model, complex, params_nested)
      }
      else{
        if (!is.null(params_master$a)) params_nested$alpha <- params_master$a else params_nested$alpha = -1
        result <- lm.logpost.bas2(y, x, model, complex, params_nested)
      }

    } else {
      stop("Unexpected error in prior_beta logic for Gaussian.")
    }
  } else {
    stop("Unsupported family: ", params_master$family, ". Supported families are 'binomial', 'poisson', 'gamma', or 'gaussian'.")
  }
  
  # Ensure consistent return structure
  if (is.null(result$crit) || is.null(result$coefs)) {
    stop("Error in computation: Returned result does not contain 'crit' and 'coefs'.")
  }
  
  return(list(crit = result$crit, coefs = result$coefs))
}