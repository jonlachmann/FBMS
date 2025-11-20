# Title     : Alpha generators
# Objective : Generation of alpha parameters for features that are kept static throughout a population
# Created by: jonlachmann
# Created on: 2021-03-16

gen.alphas <- function (strategy, feature, data, loglik, verbose) {
  if (strategy == "deep") feature <- alpha_3(feature, data, loglik, verbose)
  else if (strategy == "random") feature <- alpha_4(feature, data) else stop("Not implemented.")
  return(feature)
}

#' Alpha Generator Using Strategy 1 as per Hubin et. al.
#' TODO: This is just a placeholder.
#'
#' @noRd
#'
#' @param feature The feature to generate alphas for
alpha_1 <- function (feature) {
  return(feature)
}

#' Alpha Generator Using Strategy 2 as per Hubin et. al.
#' TODO: This is just a placeholder.
#'
#' @noRd
#'
#' @param feature The feature to generate alphas for
alpha_2 <- function (feature) {
  return(feature)
}

#' Alpha Generator Using Strategy 3 as per Hubin et. al.
#'
#' @param feature The feature to generate alphas for
#' @param data The dataset used
#' @param loglik log likelihood function to use
#'
#' @noRd
#'
alpha_3 <- function (feature, data, loglik, verbose) {
  # Create the string representation of the feature with variable alphas
  featfun <- print.feature(feature, dataset = TRUE, fixed = data$fixed, alphas = TRUE)
  featfun <- set_alphas(featfun)
  # Return if there are no alphas to set
  if (featfun$count == 0) return(feature)

  # Set initial range for Simulated Annealing
  range <- 10
  done <- FALSE
  while (!done) {
    # Run simulated annealing on current range
    sares <- GenSA::GenSA(rnorm(featfun$count), loglik,
                      rep(-range / 2, featfun$count), rep(range / 2, featfun$count),
                      control = list(max.call = 5e3, maxit = 20,nb.stop.improvement = 5), data, featfun$formula)
    # Check if any estimate is on the edge of the range, if so, extend the range and run again
    if (sum((sares$par == (-range / 2)) + (sares$par == (range / 2))) != 0) range <- range*2
    else done <- TRUE
  }
  if (sum(sares$par == 0) == featfun$count) {
    if (verbose) cat("All zero feature occured.\n")
    return(NULL)
  }
  # Inject the new alphas into the feature
  feature <- update.alphas(feature, sares$par)
  return(feature)
}

#' Alpha Generator Using Strategy 4 as per Hubin et. al.
#'
#' @param feature The feature to generate alphas for
#'
#' @noRd
#'
alpha_4 <- function (feature, data) {
  # Create the string representation of the feature with variable alphas
  featfun <- print.feature(feature, dataset = TRUE, fixed = data$fixed, alphas = TRUE)
  featfun <- set_alphas(featfun)
  
  feature <- update.alphas(feature, rnorm(featfun$count,0,1))
  return(feature)
}
