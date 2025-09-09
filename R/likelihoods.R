# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24


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
#' logistic.loglik(as.integer(rnorm(100) > 0), matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export logistic.loglik
logistic.loglik <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5))) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r = 1 / dim(x)[1])
  else if(length(mlpost_params$r) == 0)
    mlpost_params$r = 1 / dim(x)[1]
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = binomial())})
  
  if (length(mod) == 0 || is.nan(mod$deviance)) {
    return(list(crit = -.Machine$double.xmax + log_prior(mlpost_params, complex), coefs = rep(0, sum(model))))
  }
  
  
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(mlpost_params$r) * sum(complex$oc))) / 2
  return(list(crit = ret, coefs = mod$coefficients))
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
#' @importFrom stats Gamma poisson
#'
#' @examples
#' glm.loglik(abs(rnorm(100))+1, matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export glm.loglik
glm.loglik <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5), family = "Gamma")) {
  if (length(mlpost_params) == 0)
    mlpost_params <- list(r = 1 / dim(x)[1])
  else if(length(mlpost_params$r) == 0)
    mlpost_params$r = 1 / dim(x)[1]

  if (mlpost_params$family == "binomial") {
    fam = binomial()
  } else if (mlpost_params$family == "poisson") {
    fam = poisson()
  } else {
    fam = Gamma()
  }

  #browser()
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = fam)})
  
  if (length(mod) == 0 || is.nan(mod$deviance)) {
    return(list(crit = -.Machine$double.xmax + log_prior(mlpost_params, complex), coefs = rep(0, sum(model))))
  }
  
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(mlpost_params$r) * sum(complex$oc))) / 2
  return(list(crit = ret, coefs = mod$coefficients))
}

#' Log likelihood function for glm regression with Zellner's g-prior and BIC-like approximations
#'
#' This function estimates marginal likelihood for generalized linear models using a
#' BIC-style penalty adjusted to approximate Zellner's g-prior effect.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model A logical vector indicating which features are included in the model
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, including:
#'   \itemize{
#'     \item \code{r} - scalar tuning parameter for the prior (default is 1 / number of rows of \code{x})
#'     \item \code{family} - GLM family as string ("binomial", "poisson", "Gamma"), default is "binomial"
#'     \item \code{g} - scalar specifying the g prior hyperparameter (default max of model size squared and sample size)
#'   }
#' @importFrom stats Gamma poisson
#'
#' @return A list with the approximate log marginal likelihood (\code{crit}) and the posterior mode of coefficients (\code{coefs})
#'
#' @examples
#' glm.loglik.g(as.integer(rnorm(100) > 0), 
#' cbind(1, matrix(rnorm(100))), c(TRUE, TRUE), list(oc = 1),
#'  list(r = 1/100, family = "binomial", g = 10))
#'
#' @export
glm.loglik.g <- function(y, x, model, complex, mlpost_params = list(r = NULL, family = "binomial", g = NULL)) {
  if (sum(model) == 0) return(list(crit = -Inf, coefs = numeric()))
  n <- nrow(x)
  if (is.null(mlpost_params$r)) mlpost_params$r <- 1 / n
  if (is.null(mlpost_params$g)) mlpost_params$g <- max(sum(model)^2, n)
  
  # Safely handle family argument with fallback default
  family_use_str <- if (is.null(mlpost_params$family) || length(mlpost_params$family) != 1 || is.na(mlpost_params$family)) {
    "binomial"
  } else {
    mlpost_params$family
  }
  
  family_use <- switch(family_use_str,
                       binomial = stats::binomial(),
                       poisson = stats::poisson(),
                       Gamma = stats::Gamma(),
                       stats::binomial())
  
  Xsub <- as.matrix(x[, model, drop = FALSE])
  
  fit <- tryCatch({
    suppressWarnings(stats::glm.fit(Xsub, y, family = family_use))
  }, error = function(e) {
    warning("glm.loglik.g: GLM fit failed: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(fit)) return(list(crit = -Inf, coefs = numeric()))
  
  dev <- fit$deviance
  p <- fit$rank
  
  penalty_g <- p * log(1 + mlpost_params$g)
  
  crit <- (-1 / 2) * (dev + p * log(n) - 2 * log(mlpost_params$r) * sum(complex$oc) + penalty_g)
  
  list(crit = crit, coefs = fit$coefficients)
}


#' Log likelihood function for glm regression with parameter priors from BAS package
#'
#' This is a placeholder version of the function.
#' It falls back to \code{glm.loglik.g} and raises a warning if the full function is not loaded.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model A logical vector indicating which features are included in the model
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the approximate log marginal likelihood combined with the log prior (\code{crit}) and the posterior mode of coefficients (\code{coefs})
#'
#' @examples
#' glm.logpost.bas(as.integer(rnorm(100) > 0), 
#' cbind(1, matrix(rnorm(100))), c(TRUE, TRUE), 
#' list(oc = 1))
#'
#' @export
glm.logpost.bas <- function(y, x, model, complex, mlpost_params = NULL) {
  #warning("Full glm.logpost.bas not loaded; using glm.loglik.g fallback approximation")
  glm.loglik.g(y, x, model, complex, mlpost_params)
}


#' Log likelihood function for Gaussian regression with parameter priors from BAS package
#'
#' This is a placeholder version of the function.
#' It falls back to \code{gaussian.loglik.g} and raises a warning if the full function is not loaded.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model A logical vector indicating which features are included in the model
#' @param complex A list of complexity measures for the features
#' @param mlpost_params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the approximate log marginal likelihood combined with the log prior (\code{crit}) and the posterior mode of coefficients (\code{coefs})
#'
#' @examples
#' lm.logpost.bas(rnorm(100), cbind(1, matrix(rnorm(100))), c(TRUE, TRUE), list(oc = 1))
#'
#' @export
lm.logpost.bas <- function(y, x, model, complex, mlpost_params = NULL) {
  #warning("Full lm.logpost.bas not loaded; using gaussian.loglik.g fallback approximation")
  gaussian.loglik.g(y, x, model, complex, mlpost_params)
}

#' @keywords internal
#' @noRd
loadGlmBas <- function() {
  glm_code <- "
  function (y, x, model, complex, mlpost_params = list(r = NULL, family = 'binomial', beta_prior = Jeffreys(), laplace = FALSE)) {
    if (length(mlpost_params) == 0)
      mlpost_params <- list(r = 1 / dim(x)[1], family = 'binomial', beta_prior = g.prior(max(dim(x)[1], sum(model) - 1)), laplace = FALSE)
    else if(length(mlpost_params$r) == 0)
      mlpost_params$r = 1 / dim(x)[1]
    if(length(mlpost_params$laplace) == 0)
      mlpost_params$laplace = FALSE
    p <- sum(model) - 1 
    if (p == 0) {
      probinit <- as.numeric(c(1, 0.99))
      model[2] <- TRUE
    } else {
      probinit <- as.numeric(c(1, rep(0.99, p)))
    }
    mod <- NULL
    if (mlpost_params$family == 'binomial')
      family_use <- binomial()
    else if (mlpost_params$family == 'poisson')
      family_use <- poisson()
    else
      family_use <- Gamma()
    tryCatch({ suppressWarnings({
        mod <- .Call(
          BAS:::C_glm_deterministic,
          y = as.numeric(y),
          X = as.matrix(x[, model]),
          Roffset = as.numeric(rep(0, length(y))),
          Rweights = as.numeric(rep(1, length(y))),
          Rprobinit = probinit,
          Rmodeldim = as.integer(rep(0, ifelse(p == 0,2,1))),
          modelprior = uniform(),
          betaprior = mlpost_params$beta_prior,
          family = family_use,
          Rcontrol = glm.control(),
          Rlaplace =  as.integer(mlpost_params$laplace)
        )
      })
    }, error = function(e) {
      mod <- NULL
      cat('An error occurred:', conditionMessage(e), '\\n')
    })
    if (length(mod) == 0 || is.nan(mod$logmarg[2])) {
      return(list(crit = -.Machine$double.xmax + log_prior(mlpost_params, complex), coefs = rep(0, p + 1)))
    }
    if (p == 0) {
      ret <- mod$logmarg[2] + log(mlpost_params$r) * sum(complex$oc)
      return(list(crit = ret, coefs = mod$mle[[2]]))
    }
    ret <- mod$logmarg + log(mlpost_params$r) * sum(complex$oc)
    return(list(crit = ret, coefs = mod$mle[[1]]))
  }
  "
  eval(parse(text = glm_code))
}

#' @keywords internal
#' @noRd
loadLmBas <- function() {
  lm_code <- "
  function (y, x, model, complex, mlpost_params = list(r = exp(-0.5), beta_prior = list(method = 1))) {
    if (length(mlpost_params) == 0)
      mlpost_params <- list(
        r = 1 / dim(x)[1],
        beta_prior = list(method = 0, alpha = max(dim(x)[1], sum(model)^2))
      ) else if(length(mlpost_params$r) == 0) mlpost_params$r = 1 / dim(x)[1]
    p <- sum(model) - 1
    if (p == 0) {
      probinit <- as.numeric(c(1, 0.99))
      model[2] <- TRUE
    } else {
      probinit <- as.numeric(c(1, rep(0.99, p)))
    }
    mod <- NULL
    tryCatch({
        suppressWarnings({
          mod <- .Call(BAS:::C_deterministic,
                       y = y, X = as.matrix(x[, model]),
                       as.numeric(rep(1, length(y))),
                       probinit,
                       as.integer(rep(0, ifelse(p == 0,2,1))),
                       incint = as.integer(FALSE),
                       alpha = ifelse(length(mlpost_params$beta_prior$a) > 0, as.numeric(mlpost_params$beta_prior$a), NULL),
                       method = as.integer(mlpost_params$beta_prior$method),
                       modelprior = uniform(),
                       Rpivot = TRUE,
                       Rtol = 1e-7)
        })
    }, error = function(e) {
      mod <- NULL
      cat('An error occurred:', conditionMessage(e), '\\n')
    })
    if (length(mod) == 0 || is.nan(mod$logmarg[2])) {
      return(list(crit = -.Machine$double.xmax + log_prior(mlpost_params, complex), coefs = rep(0, p + 1)))
    }
    if (p == 0) {
      ret <- mod$logmarg[2] + log(mlpost_params$r) * sum(complex$oc)
      return(list(crit = ret, coefs = mod$mle[[2]]))
    }
    ret <- mod$logmarg + log(mlpost_params$r) * sum(complex$oc)
    return(list(crit = ret, coefs = mod$mle[[1]]))
  }
  "
  eval(parse(text = lm_code))
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
gaussian.loglik <- function (y, x, model, complex, mlpost_params) {
  if (sum(model) == 0)
    return(list(crit = -Inf, coefs = numeric()))
  if (length(mlpost_params) == 0)
    mlpost_params <- list()
  if (length(mlpost_params$r) == 0)
    mlpost_params$r <- 1/dim(x)[1]
  if (length(mlpost_params$beta_prior$var) == 0)
    mlpost_params$beta_prior$var <- "unknown"
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())})

  if (mlpost_params$beta_prior$var == "unknown")
    ret <- (-(mod$aic + (log(length(y)) - 2) * (mod$rank) - 2 * log(mlpost_params$r) * (sum(complex$oc)))) / 2
  else
    ret <- (-(mod$deviance / mlpost_params$beta_prior$var + log(length(y)) * (mod$rank - 1) - 2 * log_prior(mlpost_params, complex))) / 2

  return(list(crit = ret, coefs = mod$coefficients))
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
#' gaussian.loglik.g(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @export gaussian.loglik.g
gaussian.loglik.g <- function (y, x, model, complex, mlpost_params = NULL) {
  if (sum(model) == 0)
    return(list(crit = -Inf, coefs = numeric()))
  if (length(mlpost_params) == 0)
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

  if (length(mlpost_params$r) == 0 || length(mlpost_params$g) == 0) {
    mlpost_params$r <- 1 / dim(x)[1]
    if (!is.null(mlpost_params$beta_prior$g))
      mlpost_params$g <- mlpost_params$beta_prior$g
    else
      mlpost_params$g <- max(mod$rank^2, length(y))
  }

  # logarithm of marginal likelihood
  mloglik <- 0.5 * (log(1.0 + mlpost_params$g) * (dim(x)[1] - mod$rank) - log(1.0 + mlpost_params$g * (1.0 - Rsquare)) * (dim(x)[1]  - 1)) * (mod$rank != 1)

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
#' gaussian_tcch_log_likelihood(rnorm(100), matrix(rnorm(100)), c(TRUE), list(oc=1))
#'
#' @importFrom BAS phi1 hypergeometric1F1 hypergeometric2F1
#' @importFrom tolerance F1
#' @export
gaussian_tcch_log_likelihood <- function(y, x, model, complex, mlpost_params = list(r = exp(-0.5), beta_prior = list(type = "intrinsic"))) {
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
  hyper <- mlpost_params$beta_prior
  if (mlpost_params$beta_prior$type == "CH") {
    # CH prior: b and s should be user-specified, with defaults if not provided
    a <- ifelse(!is.null(hyper$a), hyper$a, 1)  # Default to 1 if not specified
    b <- ifelse(!is.null(hyper$b), hyper$b, 2)  # Default to 1 if not specified
    r <- 0
    s <- ifelse(!is.null(hyper$s), hyper$s, 1)  # Default to 1 if not specified
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "hyper-g") {
    a <- 1
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "uniform") {
    a <- 2
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "Jeffreys") {
    a <- 0.0001
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "beta.prime") {
    a <- 1/2
    b <- n - p_M - 1.5
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "benchmark") {
    a <- 0.02
    b <- 0.02 * max(n, p_M^2)
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "TG") {
    a <- 2 * ifelse(!is.null(hyper$a), hyper$a, 1)
    b <- 2
    r <- 0
    s <- 2 * ifelse(!is.null(hyper$s), hyper$s, 1)
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "ZS-adapted") {
    a <- 1
    b <- 2
    r <- 0
    s <- n + 3
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "robust") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- (n + 1) / (p_M + 1)
    k <- 1
  } else if (mlpost_params$beta_prior$type == "hyper-g-n") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- 1
    k <- 1
  } else if (mlpost_params$beta_prior$type == "intrinsic") {
    a <- 1
    b <- 1
    r <- 1
    s <- 0
    v <- (n + p_M + 1) / (p_M + 1)
    k <- (n + p_M + 1) / n
  } else if (mlpost_params$beta_prior$type == "tCCH") {
    a <- hyper$a
    b <- hyper$b
    r <- hyper$rho
    s <- hyper$s
    v <- hyper$v
    k <- hyper$k
  } else {
    stop("Unknown prior name: ", mlpost_params$beta_prior$type)
  }

  if (!is.null(r) & r == 0) {
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- phi1(b / 2, (n - 1) / 2, (a + b + p_M) / 2, s / (2 * v), min(0.8, R2_M / (v - (v - 1) * R2_M), log = TRUE))

    if (R2_M / (v - (v - 1) * R2_M) > 0.8) {
      warning("Infinite marginal log likelihood! phi1 last argument reduced to 0.8. Use a different prior_beta (Robust, Hyper-g/n, Intrinsic, or g-prior)")
    }

    term3 <- lbeta(a / 2, b / 2)
    term4 <- hypergeometric1F1(b / 2, (a + b) / 2, s / (2 * v), log = TRUE)
    marginal_likelihood <- log_likelihood + (term1) + (term2) - (p_M / 2) * log(v) - ((n - 1) / 2) * log(1 - (1 - 1 / v) * R2_M) - (term3) - (term4)
  } else if (!is.null(s) & s == 0) {
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- hypergeometric2F1(r, b / 2, (a + b) / 2, 1 - k, log = TRUE)
    term3 <- F1((a + p_M) / 2, (a + b + p_M + 1 - n - 2 * r) / 2, (n - 1) / 2, (a + b + p_M) / 2, 1 - k, 1 - k - (R2_M^2 * k) / ((1 - R2_M) * v))
    marginal_likelihood <- log_likelihood + (a + p_M - 2 * r) / 2 * log(k) + (term1) - (term2) - (term3) - (p_M / 2) * log(v) - log(1 - R2_M) * ((n - 1) / 2) - lbeta(a / 2, b / 2)

  } else {
    stop("Invalid inputs: either r = 0 or s = 0 must be specified.")
  }

  if (length(mlpost_params$r) == 0) mlpost_params$r <- 1 / dim(x)[1]  # default value or parameter r

  lp <- log_prior(mlpost_params, complex)

  if (is.nan(marginal_likelihood)) {
    return(list(crit = -.Machine$double.xmax + lp, coefs = rep(0,sum(model))))
  }
  
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
#' logistic.loglik.ala(as.integer(rnorm(100) > 0), matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export logistic.loglik.ala
logistic.loglik.ala <- function (y, x, model, complex, mlpost_params = list(r = exp(-0.5))) {
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
#' @export logistic.loglik.alpha
logistic.loglik.alpha <- function (a, data, mu_func) {
  m <- 1 / (1 + exp(-eval(parse(text = mu_func))))
  -sum((data$y[,1] * log(m) + (1 - data$y[, 1]) * log(1 - m)))
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
#'gaussian.loglik.alpha(my_alpha,my_data,my_mu)
#'}
#' @export gaussian.loglik.alpha
gaussian.loglik.alpha <- function (a, data, mu_func) {
  m <- eval(parse(text = mu_func))
  sum((data$y[, 1] - m)^2)
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
#'   - a: Tuning parameter for g-prior (default: max(n, p^2))
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
#' fbms.mlik.master(y = rnorm(100), 
#' x = matrix(rnorm(100)), 
#' c(TRUE,TRUE), 
#' list(oc = 1),
#' mlpost_params = list(family = "gaussian", beta_prior = list(type = "g-prior", a = 2),
#'          r = exp(-0.5)))
#'
#' @importFrom BAS robust beta.prime bic.prior CCH EB.local g.prior hyper.g hyper.g.n tCCH intrinsic TG Jeffreys uniform
#' @export
fbms.mlik.master <- function(y, x, model, complex, mlpost_params = list(family = "gaussian", beta_prior = list(type = "g-prior"), r = NULL)) {
  # Extract dimensions
  n <- length(y)
  p <- length(model) - 1  # Number of predictors excluding intercept
  params_use <- list()

  if(length(mlpost_params$r) == 0)
    mlpost_params$r <-  1/length(y)
  
  if(mlpost_params$family == "gaussian")
    params_use$beta_prior <- gen.mlpost.params.lm(mlpost_params$beta_prior$type, mlpost_params$beta_prior, p, n)
  else
  {
    params_use$beta_prior <- gen.mlpost.params.glm(mlpost_params$beta_prior$type, mlpost_params$beta_prior, p, n)
    params_use$family <- mlpost_params$family
  }
  
  params_use$r <- mlpost_params$r
  
  
  loglik.pi <- select.mlpost.fun(mlpost_params$beta_prior$type, mlpost_params$family)
  
  
  result <- loglik.pi(y, x, model, complex, params_use)

  return(list(crit = result$crit, coefs = result$coefs))
}