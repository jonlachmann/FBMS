# Title     : Alpha generators
# Objective : Generation of alpha parameters for features that are kept static throughout a population
# Created by: jonlachmann
# Created on: 2021-03-16

gen.alphas <- function (strategy, feature, data, loglik) {
  if (strategy == 1) stop("Not implemented.")
  else if (strategy == 2) stop("Not implemented.")
  else if (strategy == 3) feature <- alpha_3(feature, data, loglik)
  return(feature)
}

#' Alpha generator using strategy 1 as per Hubin et. al.
#' TODO: Implement me
#' @param feature The feature to generate alphas for
#' @param transforms The transforms used
alpha_1 <- function (feature) {
  return(feature)
}

#' Alpha generator using strategy 2 as per Hubin et. al.
#' TODO: Implement me
#' @param feature The feature to generate alphas for
#' @param transforms The transforms used
alpha_2 <- function (feature) {
  return(feature)
}

#' Alpha generator using strategy 3 as per Hubin et. al.
#'
#' @param feature The feature to generate alphas for
#' @param data The dataset used
#' @param transforms The transforms used
#' @param loglik log likelihood function to use
alpha_3 <- function (feature, data, loglik) {
  # Create the string representation of the feature with variable alphas
  featfun <- print.feature(feature, dataset = T, alphas = T)
  featfun <- set_alphas(featfun)
  # Return if there are no alphas to set
  if(featfun$count == 0) return(feature)

  # Set initial range for Simulated Annealing
  cat("Generating alphas\n")
  range <- 10
  done <- FALSE
  while(!done) {
    # Run simulated annealing on current range
    sares <- GenSA(rnorm(featfun$count), loglik,
                      rep(-range/2,featfun$count), rep(range/2,featfun$count),
                      control=list(max.call=5e3), data, featfun$formula)
    # Check if any estimate is on the edge of the range, if so, extend the range and run again
    if (sum((sares$par==(-range/2))+(sares$par==(range/2))) != 0) range <- range*2
    else done <- TRUE
  }
  if (sum(sares$par==0) == featfun$count) {
    cat("All zero feature occured.\n")
    return(NULL)
  }
  # Inject the new alphas into the feature
  feature <- update.alphas(feature, sares$par)
  return(feature)
}