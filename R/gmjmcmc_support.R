# Title     : GMJMCMC Support functions
# Objective : Support functions for GMJMCMC (Genetically Modified MJMCMC) algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

#' Set the transformations option for GMJMCMC (Genetically Modified MJMCMC),
#' this is also done when running the algorithm, but this function allows for it to be done manually.
#'
#' @param transforms The vector of non-linear transformations
#'
#' @return No return value, just sets the gmjmcmc-transformations option
#'
#' @examples
#' set.transforms(c("p0","p1"))
#' 
#'
#' @export set.transforms
set.transforms <- function (transforms) {
  old_transforms <- getOption("gmjmcmc-transformations")
  options("gmjmcmc-transformations" = transforms)
  return(old_transforms)
}

# Function to verify inputs and help the user find if they did anything wrong
# TODO: Finish this!!
verify.inputs <- function (data, loglik.pi, transforms, T, N, N.final, probs, params) {
  # Get information about the data
  n.cov <- ncol(data) - 1
  n.obs <- nrow(data)

  # Get information about the transforms
  n.trans <- length(transforms)

  # Verify that the probability list is formatted as it should be
  if (length(probs$large) != 1) error <- c(error, "Too many large jump probabilities (should be 1).")
  else if (probs$large > 1 | probs$large < 0) error <- c(error, "Large jump probability must be in [0,1].")
  if (length(probs$large.kern) != 4) error <- c(error, "There should be 4 large jump kernel probabilities.")
  if (sum(probs$large.kern > 1 | probs$large.kern < 0) != 0) error <- c(error, "Large jump kernel probabilities must be in [0,1].")
}

#' Function for calculating marginal inclusion probabilities of features given a list of models
#' @param models The list of models to use.
#'
#' @return A numeric vector of marginal model probabilities based on relative frequencies of model visits in MCMC.
#'
#' @examples
#' result <- gmjmcmc(matrix(rnorm(600), 100), P = 2, gaussian.loglik, NULL, c("p0", "exp_dbl"))
#' marginal.probs(result$models[[1]])
#'
#' @export
marginal.probs <- function (models) {
  mod.count <- length(models)
  probs <- rep(0, length = length(models[[1]]$model))
  for (i in 1:mod.count) {
    probs <- probs + models[[i]]$model
  }
  probs <- probs / mod.count
  return(probs)
}

#' Function for calculating feature importance through renormalized model estimates
#' @param models The models to use.
#' @param type Select which probabilities are of interest, features or models
#' 
#' @noRd
marginal.probs.renorm <- function (models, type = "features") {
  models <- lapply(models, function (x) x[c("model", "crit")])
  model.size <- length(models[[1]]$model)
  models.matrix <- matrix(unlist(models), ncol = model.size + 1, byrow = TRUE)
  duplicates <- duplicated(models.matrix[, 1:(model.size)], dim = 1, fromLast = TRUE)
  models.matrix <- models.matrix[!duplicates, ]
  if(!is.matrix(models.matrix))
    models.matrix <- t(as.matrix(models.matrix))
  
  max_mlik <- max(models.matrix[,(model.size + 1)])
  crit.sum <- sum(exp(models.matrix[, (model.size + 1)] - max_mlik))
  if (type == "features" || type == "both") {
    probs.f <- matrix(NA,1, model.size)
    for (i in 1:(model.size)) probs.f[i] <- sum(exp(models.matrix[as.logical(models.matrix[, i]),(model.size + 1)] - max_mlik)) / crit.sum
  }
  if (type == "models" || type == "both") {
    probs.m <- matrix(NA,1, nrow(models.matrix))
    for (i in seq_len(nrow(models.matrix))) probs.m[i] <- exp(models.matrix[i, ncol(models.matrix)] - max_mlik) / crit.sum
  }

  if (type == "features") {
    result <- list(idx = which(!duplicates), probs = probs.f)
  } else if (type == "models") {
    result <- list(idx = which(!duplicates), probs = probs.m)
  } else {
    result <- list(idx = which(!duplicates), probs.f = probs.f, probs.m = probs.m)
  }
  return(result)
}


# Function for precalculating features for a new feature population
precalc.features <- function (data, features) {
  precalc <- matrix(NA, nrow(data), length(features) + 2)
  precalc[, 1:2] <- data[, 1:2]
  for (f in seq_along(features)) {
    feature_string <- print.feature(features[[f]], dataset = TRUE)
    precalc[, (f + 2)] <- eval(parse(text = feature_string))
  }
  # Replace any -Inf and Inf values caused by under- or overflow
  precalc <- replace.infinite.data.frame(precalc)
  return(precalc)
}

# TODO: Compare to previous mliks here instead, also add a flag to do that in full likelihood estimation scenarios.
# Function to call the model function
loglik.pre <- function (loglik.pi, model, complex, data, params = NULL, visited.models = visited.models, sub = sub) {
  if (!is.null(visited.models) && has_key(visited.models, model)) {
    if (!sub) {
      return(visited.models[[model]])
    } else {
      params$coefs <- visited.models[[model]]$coefs
      params$crit <- visited.models[[model]]$crit
    }
  }
  # Get the complexity measures for just this model
  complex <- list(width = complex$width[model], oc = complex$oc[model], depth = complex$depth[model])
  # Call the model estimator with the data and the model, note that we add the intercept to every model
  model.res <- loglik.pi(data[, 1], data[, -1], c(T, model), complex, params)
  # Check that the critical value is acceptable
  if (!is.numeric(model.res$crit) || is.nan(model.res$crit)) model.res$crit <- -.Machine$double.xmax
  # Alpha cannot be calculated if the current and proposed models have crit which are -Inf or Inf
  if (is.infinite(model.res$crit)) {
    if (model.res$crit > 0) model.res$crit <- .Machine$double.xmax
    else model.res$crit <- -.Machine$double.xmax
  }
  return(model.res)
}

# Function to check the data
# Checks that there is an intercept in the data, adds it if missing
# Coerces the data to be of type matrix
check.data <- function (data, verbose) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
    if (verbose) cat("Data coerced to matrix type.\n")
  }
  if (sum(data[, 2] == 1) != nrow(data)) {
    data <- cbind(data[, 1], 1, data[, -1])
    if (verbose) cat("Intercept added to data.\n")
  }
  return(data)
}

# Function to get the dimensions of a dataset, adding an intercept if necessary
data.dims <- function (data) {
  dims <- dim(data)
  if (sum(data[,2] == 1) != nrow(data)) {
    dims[2] <- dims[2] + 1
  }
  return(dims)
}

# Function to extract column names if they are well formed
get.labels <- function (data, verbose) {
  labels <- colnames(data)[-(1:2)]
  if (is.null(labels)) return(F)
  if (sum(is.na(labels)) != 0) {
    if (verbose) cat("NA labels present, using x#\n")
    return(F)
  }
  return(labels)
}