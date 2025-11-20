#' Fit a BGNLM Model Using MJMCMC or GMJMCMC Sampling.
#'
#' This function fits a model using the relevant MCMC sampling. The user can specify the formula,
#' family, data, transforms, and other parameters to customize the model.
#'
#' @param formula A formula object specifying the model structure. Default is NULL.
#' @param family The distribution family of the response variable. Currently supports "gaussian", "binomial", "poisson", "gamma", and  "custom". Default is "gaussian".
#' @param beta_prior Type of prior as a string (default: "g-prior" with a = max(n, p^2)). Possible values include:
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
#' @param model_prior a list with parameters of model priors, by default r should be provided
#' @param extra_params extra parameters to be passed to the loglik.pi function
#' @param loglik.pi Custom function to compute the logarithm of the posterior mode based on logarithm of marginal likelihood and logarithm of prior functions (needs specification only used if family = "custom")
#' @param data A data frame or matrix containing the data to be used for model fitting. If the outcome variable is in the first column of the data frame, the formula argument in fbms can be omitted, provided that all other columns are intended to serve as input covariates.
#' @param method Which fitting algorithm should be used, currently implemented options include "gmjmcmc", "gmjmcmc.parallel", "mjmcmc" and "mjmcmc.parallel" with "mjmcmc" being the default and 'mjmcmc' means that only linear models will be estimated
#' @param verbose If TRUE, print detailed progress information during the fitting process. Default is TRUE.
#' @param impute TRUE  means imputation combined with adding a dummy column with indicators of imputed values, FALSE (default) means only full data is used.
#' @param ... Additional parameters to be passed to the underlying method.
#'
#' @return An object containing the results of the fitted model and MCMC sampling.
#'
#' @examples
#' # Fit a Gaussian multivariate time series model
#' fbms_result <- fbms(
#'  X1 ~ .,
#'  family = "gaussian",
#'  method = "gmjmcmc.parallel",
#'  data = data.frame(matrix(rnorm(600), 100)),
#'  transforms = c("sin","cos"),
#'  P = 10,
#'  runs = 1,
#'  cores = 1
#' )
#' summary(fbms_result)
#' 
#'
#' @seealso \code{\link{mjmcmc}}, \code{\link{gmjmcmc}}, \code{\link{gmjmcmc.parallel}}
#' @export
#' @importFrom stats terms
fbms <- function (
  formula = NULL,
  family = "gaussian",
  beta_prior = list(type = "g-prior"),
  model_prior = NULL,
  extra_params = NULL,
  data = NULL,
  impute = FALSE,
  loglik.pi = NULL,
  method = "mjmcmc",
  verbose = TRUE,
  ...
) {
  if (length(data) == 0)
    stop("Training data must be provided!")

  if (length(model_prior) == 0)
    model_prior = list(r = 1 / dim(data)[1])
  if (family != "custom") {
    mlpost_params <- model_prior
    loglik.pi <- select.mlpost.fun(beta_prior$type, family)
    if (family == "gaussian") {
      mlpost_params$beta_prior <- gen.mlpost.params.lm(beta_prior$type, beta_prior, ncol(data) - 1, nrow(data))
    } else {
      mlpost_params$beta_prior <- gen.mlpost.params.glm(beta_prior$type, beta_prior, ncol(data) - 1, nrow(data))
      mlpost_params$beta_prior$type <- beta_prior$type
      mlpost_params$family <- family
    }
  } else {
    loglik.pi <- loglik.pi
    mlpost_params <- c(model_prior, beta_prior, extra_params)
  }

  if (!is.null(formula)) {
    if (missing(data)) {
      data <- environment(formula)
    }

    na.opt <- getOption("na.action")
    if (impute)
      options(na.action = 'na.pass')
    else
      options(na.action = 'na.omit')
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    Y <- model.response(mf, "any")
    X <- model.matrix(formula, data = data)
    intercept <- attr(terms(formula, data = data), "intercept") == 1
    if (intercept) X <- X[, -1, drop = FALSE]
    mis.Y <- which(is.na(Y))
    if (length(mis.Y) > 0) {
      warning("Missing values in the response. Dropped.")
      Y <- Y[-c(mis.Y)]
      X <- X[-c(mis.Y), ]
    }
    
    mis.X <- sum(is.na(X))
    imputed <- NULL
    if (impute && mis.X > 0) {
      print("Imputing missing values!")
      X <- data.frame(X)
      na.matr <- data.frame(1 * (is.na(X)))
      names(na.matr) <- paste0("mis_", names(na.matr))
      cm <- colMeans(na.matr)
      na.matr <- na.matr[, cm != 0]
      for (i in seq_along(X)) {
        X[[i]][is.na(X[[i]])] <- median(X[[i]], na.rm = TRUE)
      }
      imputed <- names(X)[cm != 0]
      X <- data.frame(X, na.matr)
      rm(na.matr)
      rm(cm)
      print("Continue to sampling!")
    } else if (mis.X > 0) {
      print("Dropping missing values!")
    }
  } else {
    Y <- data[, 1]
    X <- data[, -1, drop = FALSE]
    intercept <- TRUE
    imputed <- NULL
    na.opt <- getOption("na.action")
    if (impute) {
      options(na.action = 'na.pass')
      stop("Imputation is only implemented when formula is provided.\n Please specify formula and rerun!")
    } else {
      options(na.action = 'na.omit')
    }
  }
  
  if (method == "mjmcmc.parallel")
    res <- mjmcmc.parallel(x = X, y = Y, loglik.pi = loglik.pi, mlpost_params = mlpost_params, intercept = intercept, verbose = verbose, ...)
  else if (method == "mjmcmc")
    res <- mjmcmc(x = X, y = Y, loglik.pi = loglik.pi, mlpost_params = mlpost_params, intercept = intercept, verbose = verbose, ...)
  else if (method == "gmjmcmc.parallel")
    res <- gmjmcmc.parallel(x = X, y = Y, loglik.pi = loglik.pi, mlpost_params = mlpost_params, intercept = intercept, verbose = verbose,...)
  else if (method == "gmjmcmc")
    res <- gmjmcmc(x = X, y = Y, loglik.pi = loglik.pi, mlpost_params = mlpost_params, intercept = intercept, verbose = verbose, ...)
  else
    stop("Error: Method must be one of gmjmcmc, gmjmcmc.parallel, mjmcmc or mjmcmc.parallel!")
  
  attr(res, "imputed") <- imputed
  attr(res, "all_names") <- names(X)[1:(dim(X)[2] - 1)]
  options(na.action = na.opt)
  return(res)
}

gen.mlpost.params.glm <- function (beta_prior, user_params, p, n) {
  
  if(beta_prior == "Jeffreys-BIC")
  {    
    return(NULL)
  } else if(beta_prior == "beta.prime") {
    return(BAS::beta.prime(n = n))
  } else if (beta_prior == "CH") {
    check_required_params(c("a", "b", "s"), user_params, beta_prior)
    return(BAS::CCH(alpha = user_params$a, beta = user_params$b, s = user_params$s))
  } else if (beta_prior == "EB-local") {
    return(BAS::EB.local())
  } else if (beta_prior == "g-prior") {
    if (is.null(user_params$g)) {
      user_params$g <- max(p^2, n)
    }
    return(BAS::g.prior(user_params$g))
  } else if (beta_prior == "hyper-g") {
    check_required_params("a", user_params, beta_prior)
    params <- BAS::hyper.g(alpha = user_params$a)
    params$method <- 1
    return(params)
  } else if (beta_prior == "tCCH") {
    check_required_params(c("a", "b", "s", "rho", "v", "k"), user_params, beta_prior)
    return(BAS::tCCH(
      alpha = user_params$a,
      beta = user_params$b,
      s = user_params$s,
      r = user_params$rho,
      v = user_params$v,
      theta = user_params$k
    ))
  } else if (beta_prior == "intrinsic") {
    return(BAS::intrinsic(n = n))
  } else if (beta_prior == "TG") {
    check_required_params("a", user_params, beta_prior)
    return(BAS::TG(alpha = user_params$a))
  } else if (beta_prior == "Jeffreys") {
    return(BAS::Jeffreys())
  } else if (beta_prior == "uniform") {
    return(BAS::tCCH(alpha = 2, beta = 2, s = 0, r = 0, v = 1, theta = 1))
  } else if (beta_prior == "benchmark") {
    return(BAS::tCCH(alpha = 0.02, beta = 0.02 * max(n, p^2), s = 0, r = 0, v = 1, theta = 1))
  } else if (beta_prior == "ZS-adapted") {
    return(BAS::tCCH(alpha = 1, beta = 2, s = n + 3, r = 0, v = 1, theta = 1))
  } else if (beta_prior == "robust") {
    return(BAS::robust(n = as.numeric(n)))# important to cast to numeric for BAS, do not change.
  } else if (beta_prior == "hyper-g-n") {
    if (is.null(user_params$a)) user_params$a <- 3
    return(BAS::hyper.g.n(alpha = user_params$a, n = n))
  } else if (beta_prior == "BIC") {
    return(BAS::bic.prior(n = n))
  } 
  stop("Unknown prior, please verify your inputs.")
}


gen.mlpost.params.lm <- function (beta_prior, user_params, p, n) {
  
  if (beta_prior == "Jeffreys-BIC") {
    if(length(user_params$var)==0) 
    {  
      user_params$var <- "unknown"
    }
    return(list(var = user_params$var))
  }else if (beta_prior == "beta.prime") {
    return(list(type = "beta.prime"))
  } else if (beta_prior == "CH") {
    check_required_params(c("a", "b", "s"), user_params, beta_prior)
    user_params <- list(type =
      "CH",
      a = user_params$a,
      b = user_params$b,
      s = user_params$s
    )
    return(user_params)
  } else if (beta_prior == "tCCH") {
    check_required_params(c("a", "b", "s", "rho", "v", "k"), user_params, beta_prior)
    user_params <- list(
      type = "tCCH",
      a = user_params$a,
      b = user_params$b,
      s = user_params$s,
      rho = user_params$rho,
      v = user_params$v,
      k = user_params$k
    )
    return(user_params)
  } else if (beta_prior == "intrinsic") {
    return(list(type = "intrinsic"))
  } else if (beta_prior == "TG") {
    check_required_params(c("a", "s"), user_params, beta_prior)
    user_params <- list(
      type = "TG",
      a = user_params$a,
      s = user_params$s
    )
    return(user_params)
  } else if (beta_prior == "Jeffreys") {
    return(list(type = "Jeffreys"))
  } else if (beta_prior == "ZS-adapted") {
    return(list(type = "ZS-adapted"))
  } else if (beta_prior == "benchmark") {
    return(list(type = "benchmark"))
  } else if (beta_prior == "robust") {
    return(list(type = "robust"))
  } else if (beta_prior == "uniform") {
    return(list(type = "uniform"))
  } else{
    if(!is.null(user_params$g))
      user_params$a <- user_params$g
    if (!is.null(user_params$a)) {
      alpha <- user_params$a 
    } else {
      if (beta_prior == "g-prior") {
        alpha <- min(p^2, n)
      } else {
        alpha <- -1 #check how BAS uses the default
      }
    }
    if (beta_prior == "g-prior") {
      return(list(method = 0, alpha = alpha))
    } else if (beta_prior == "hyper-g") {
      return(list(method = 1, alpha = alpha))
    } else if (beta_prior == "EB-local") {
      return(list(method = 2, alpha = alpha))
    } else if (beta_prior == "BIC") {
      return(list(method = 3, alpha = alpha))
    } else if (beta_prior == "ZS-null") {
      return(list(method = 4, alpha = alpha))
    } else if (beta_prior == "ZS-full") {
      return(list(method = 5, alpha = alpha))
    } else if (beta_prior == "hyper-g-laplace") {
      return(list(method = 6, alpha = alpha))
    } else if (beta_prior == "AIC") {
      return(list(method = 7, alpha = alpha))
    } else if (beta_prior == "EB-global") {
      return(list(method = 2, alpha = alpha))
    } else if (beta_prior == "JZS") {
      return(list(method = 9, alpha = alpha))
    } else if (beta_prior == "hyper-g-n") {
      return(list(method = 8, alpha = alpha))
    } else {
      stop("Unrecognized prior_beta for Gaussian GLM: ", beta_prior)
    }
  }
}
  
check_required_params <- function (required, user_params, beta_prior) {
  for (req in required) {
    if (is.null(user_params[[req]]) || !is.numeric(user_params[[req]])) {
      par_names <- paste0(required, collapse = ", ")
      stop(paste0("The parameters: ", par_names, " must be provided for the ", beta_prior, " prior."))
      return(FALSE)
    }
  }
  return(TRUE)
}

select.mlpost.fun <- function (beta_prior, family) {
  if (!(family %in% c("binomial", "poisson", "gamma", "gaussian"))) {
    stop(paste0(
      "Unsupported family: ", family, ". Supported families are 'binomial', 'poisson', 'gamma', or 'gaussian'."
    ))
  }
  gaussian_only_priors <- c("ZS-null", "ZS-full", "hyper-g-laplace","BIC", "AIC", "JZS","EB-global")
  gaussian_not_robust <- c("CH", "tCCH", "ZS-adapted", "TG", "beta.prime", "benchmark", "Jeffreys")
  gaussian_robust <- c("g-prior", "hyper-g", "EB-local","BIC", "Jeffreys-BIC", "ZS-null", "ZS-full", "hyper-g-laplace",
                       "AIC",  "hyper-g-n",  "JZS")
  gaussian_tcch <- c("CH", "tCCH", "TG","beta.prime", "intrinsic", "ZS-adapted", "uniform","Jeffreys", "benchmark", "robust")
  gaussian_bas <- c("g-prior", "hyper-g", "EB-local","ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS")
  glm_priors <- c("CH", "tCCH", "TG","beta.prime", "EB-local", "g-prior", "hyper-g", "hyper-g-n",
                               "intrinsic", "ZS-adapted", "Jeffreys", "uniform", "benchmark", "robust", "Jeffreys-BIC")

  if (family %in% c("binomial", "poisson", "gamma")) {
    if (beta_prior %in% gaussian_only_priors) {
      stop(paste0(
        "Prior ", beta_prior, " is not supported for the GLM family", family,
        ". Supported priors are: ", paste(glm_priors, collapse = ", ")
      ))
    }
    if (beta_prior == "Jeffreys-BIC") {
      if (family == "binomial") {
        return(logistic.loglik)
      } else {
        return(glm.loglik)
      }
    } else {
      return(loadGlmBas())
    }
  } else if (family == "gaussian") {
    if (beta_prior %in% gaussian_not_robust) {
      warning(paste0(
        "Prior ", beta_prior, " is not recommended for Gaussian family models as it may be unstable for strong signals (R^2 > 0.9).",
        "Recommended priors under the Gaussian family are: ", paste(gaussian_robust, collapse = ", ")
      ))
    }
    if (beta_prior %in% gaussian_tcch) {
      return(gaussian_tcch_log_likelihood)
    } else if (beta_prior == "Jeffreys-BIC") {
      return(gaussian.loglik)
    } else if (beta_prior %in% gaussian_bas) {
      return(loadLmBas())
    }
  }
  stop("Unknown prior, please verify your inputs.")
}

