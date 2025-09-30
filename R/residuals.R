#' Residuals for GMJMCMC Model
#'
#' Computes residuals as the difference between observed and predicted values.
#'
#' @param object Object of class "gmjmcmc".
#' @param y Respnse.
#' @param x Covariates.
#' @param ... Additional arguments (ignored).
#' @return Vector of residuals.
#' @method residuals gmjmcmc
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "gmjmcmc", transforms = c("sigmoid"))
#' hist(residuals(model, exoplanet[,1], exoplanet[,-1]))
residuals.gmjmcmc <- function(object, y, x, ...) {
  stopifnot(inherits(object, "gmjmcmc"))
  if (is.null(object$residuals)) {
    pred <- predict(object, x)$aggr$mean
    return(y - pred)
  } else {
    object$residuals
  }
}

#' Residuals for MJMCMC Model
#'
#' Computes residuals as the difference between observed and predicted values.
#'
#' @param object Object of class "mjmcmc".
#' @param y Respnse.
#' @param x Covariates.
#' @param ... Additional arguments (ignored).
#' @return Vector of residuals.
#' @method residuals mjmcmc
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc")
#' hist(residuals(model, exoplanet[,1], exoplanet[,-1]))
residuals.mjmcmc <- function(object, y, x, ...) {
  stopifnot(inherits(object, "mjmcmc"))
  if (is.null(object$residuals)) {
    pred <- predict(object, x)$aggr$mean
    return(y - pred)
  } else {
    object$residuals
  }
}

#' Residuals for BGNLM Model
#'
#' Computes residuals as the difference between observed and predicted values.
#'
#' @param object Object of class "bgnlm_model".
#' @param y Respnse.
#' @param x Covariates.
#' @param ... Additional arguments (ignored).
#' @return Vector of residuals.
#' @method residuals bgnlm_model
#' @export
#' @examples
#' library(FBMS)
#' data(exoplanet)
#' model <- get.best.model(fbms(semimajoraxis ~ ., data = exoplanet, family = "gaussian"))
#' hist(residuals(model, exoplanet[,1], exoplanet[,-1]))
residuals.bgnlm_model <- function(object,y, x, ...) {
  stopifnot(inherits(object, "bgnlm_model"))
  if (is.null(object$residuals)) {
    pred <- predict(object, x)
    return(y - pred)
  } else {
    object$residuals
  }
}

#' Residuals for MJMCMC Parallel Model
#'
#' Computes residuals as the difference between observed and predicted values.
#'
#' @param object Object of class "mjmcmc_parallel".
#' @param y Respnse.
#' @param x Covariates.
#' @param ... Additional arguments (ignored).
#' @return Vector of residuals.
#' @method residuals mjmcmc_parallel
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc.parallel",runs = 2, cores = 1)
#' hist(residuals(model, exoplanet[,1], exoplanet[,-1]))
residuals.mjmcmc_parallel <- function(object, y, x, ...) {
  stopifnot(inherits(object, "mjmcmc_parallel"))
  if (is.null(object$residuals)) {
    pred <- predict(object, x)$aggr$mean
    return(y - pred)
  } else {
    object$residuals
  }
}


#' Residuals for GMJMCMC Merged Model
#'
#' Computes residuals as the difference between observed and predicted values.
#'
#' @param object Object of class "gmjmcmc_merged".
#' @param y Respnse.
#' @param x Covariates.
#' @param ... Additional arguments (ignored).
#' @return Vector of residuals.
#' @method residuals gmjmcmc_merged
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, 
#' method = "gmjmcmc.parallel", transforms = c("sigmoid"), 
#' runs = 2, cores = 1)
#' hist(residuals(model, exoplanet[,1], exoplanet[,-1]))
residuals.gmjmcmc_merged <- function(object, y, x, ...) {
  stopifnot(inherits(object, "gmjmcmc_merged"))
  if (is.null(object$residuals)) {
    pred <- predict(object, x)$aggr$mean
    return(y - pred)
  } else {
    object$residuals
  }
}