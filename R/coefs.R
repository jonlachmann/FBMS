#' Coefficients for GMJMCMC Model
#'
#' Extracts coefficients from the best GMJMCMC model found.
#'
#' @param object Object of class "gmjmcmc".
#' @param ... Additional arguments (ignored).
#' @return Vector of coefficients from the best model found.
#' @method coef gmjmcmc
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "gmjmcmc", transforms = c("sigmoid"))
#' coef(model)
coef.gmjmcmc <- function(object, ...) {
  stopifnot(inherits(object, "gmjmcmc"))
  cat("Posterior mode for the parameters of the best found single model:\n")
  best.mod <- get.best.model(object)
  best.mod$coefs
}

#' Coefficients for MJMCMC Model
#'
#' Extracts coefficients from the best MJMCMC model.
#'
#' @param object Object of class "mjmcmc".
#' @param ... Additional arguments (ignored).
#' @return Vector of coefficients from the best model found.
#' @method coef mjmcmc
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc")
#' coef(model)
coef.mjmcmc <- function(object, ...) {
  stopifnot(inherits(object, "mjmcmc"))
  cat("Posterior mode for the parameters of the best found single model:\n")
  best.mod <- get.best.model(object)
  best.mod$coefs
}

#' Coefficients for BGNLM Model
#'
#' Extracts coefficients from a BGNLM model.
#'
#' @param object Object of class "bgnlm_model".
#' @param ... Additional arguments (ignored).
#' @return Vector of coefficients.
#' @method coef bgnlm_model
#' @export
#' @examples
#' data(exoplanet)
#' model <- get.best.model(fbms(semimajoraxis ~ ., data = exoplanet, family = "gaussian"))
#' coef(model)
coef.bgnlm_model <- function(object, ...) {
  stopifnot(inherits(object, "bgnlm_model"))
  cat("Posterior mode for the parameters of the best found single model:\n")
  object$coefs
}

#' Coefficients for MJMCMC Parallel Model
#'
#' Extracts coefficients from the best MJMCMC parallel model.
#'
#' @param object Object of class "mjmcmc_parallel".
#' @param ... Additional arguments (ignored).
#' @return Vector of coefficients from the best model found.
#' @method coef mjmcmc_parallel
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc.parallel", cores = 1, runs = 2)
#' coef(model)
coef.mjmcmc_parallel <- function(object, ...) {
  stopifnot(inherits(object, "mjmcmc_parallel"))
  cat("Posterior mode for the parameters of the best found single model:\n")
  best.mod <- get.best.model(object)
  best.mod$coefs
}

#' Coefficients for GMJMCMC Merged Model
#'
#' Extracts coefficients from the best GMJMCMC merged model.
#'
#' @param object Object of class "gmjmcmc_merged".
#' @param ... Additional arguments (ignored).
#' @return Vector of coefficients from the best model found.
#' @method coef gmjmcmc_merged
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, 
#' method = "gmjmcmc.parallel", transforms = c("sigmoid"), 
#' runs = 2, cores = 1)
#' coef(model)
coef.gmjmcmc_merged <- function(object, ...) {
  stopifnot(inherits(object, "gmjmcmc_merged"))
  best.mod <- get.best.model(object)
  best.mod$coefs
}