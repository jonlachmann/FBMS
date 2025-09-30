
#' Print GMJMCMC Model Object
#'
#' Displays a concise summary of a GMJMCMC model object.
#'
#' @param x Object of class "gmjmcmc".
#' @param ... Additional arguments passed to summary method.
#' @return Prints a summary of the model and returns NULL
#' @method print gmjmcmc
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet,method = "gmjmcmc", transforms = c("sigmoid"))
#' print(model)
print.gmjmcmc <- function(x, ...) {
  stopifnot(inherits(x, "gmjmcmc")) 
  cat("GMJMCMC Model Summary:\n")
  print(summary(x, ...))
}

#' Print MJMCMC Model Object
#'
#' Displays a concise summary of an MJMCMC model object.
#'
#' @param x Object of class "mjmcmc".
#' @param ... Additional arguments passed to summary method.
#' @return Prints a summary of the model and returns NULL
#' @method print mjmcmc
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc")
#' print(model)
print.mjmcmc <- function(x, ...) {
  stopifnot(inherits(x, "mjmcmc"))
  cat("MJMCMC Model Summary:\n")
  print(summary(x, ...))
}

#' Print BGNLM Model Object
#'
#' Displays the coefficients of a BGNLM model object.
#'
#' @param x Object of class "bgnlm_model".
#' @param ... Additional arguments (ignored).
#' @return Prints a summary of the model and returns NULL
#' @method print bgnlm_model
#' @export
#' @examples
#' data(exoplanet)
#' model <- get.best.model(fbms(semimajoraxis ~ ., data = exoplanet, 
#' family = "gaussian"))
#' print(model)
#' model <- get.mpm.model(fbms(semimajoraxis ~ ., data = exoplanet, 
#' family = "gaussian"), y = exoplanet[,1],x = exoplanet[,-1])
#' print(model)
print.bgnlm_model <- function(x, ...) {
  stopifnot(inherits(x, "bgnlm_model"))
  cat("BGNLM Model Coefficients:\n")
  print(x$coefs)
}


#' Print MJMCMC Parallel Model Object
#'
#' Displays a concise summary of an MJMCMC parallel model object.
#'
#' @param x Object of class "mjmcmc_parallel".
#' @param ... Additional arguments passed to summary method.
#' @return Prints a summary of the model and returns NULL
#' @method print mjmcmc_parallel
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc.parallel", cores = 1, runs = 2)
#' print(model)
print.mjmcmc_parallel <- function(x, ...) {
  stopifnot(inherits(x, "mjmcmc_parallel"))
  cat("MJMCMC Parallel Model Summary:\n")
  print(summary(x, ...))
}

#' Print GMJMCMC Merged Model Object
#'
#' Displays a concise summary of a GMJMCMC merged model object.
#'
#' @param x Object of class "gmjmcmc_merged".
#' @param ... Additional arguments passed to summary method.
#' @return Prints a summary of the model and returns NULL
#' @method print gmjmcmc_merged
#' @export
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, 
#' method = "gmjmcmc.parallel", cores = 1, runs = 2, 
#' transforms = c("sigmoid"))
#' 
#' print(model)
print.gmjmcmc_merged <- function(x, ...) {
  stopifnot(inherits(x, "gmjmcmc_merged"))
  cat("GMJMCMC Merged Model Summary:\n")
  print(summary(x, ...))
}
