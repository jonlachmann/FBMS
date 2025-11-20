# Title     : GMJMCMC Support functions
# Objective : Support functions for GMJMCMC (Genetically Modified MJMCMC) algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

#' Set the Transformations Option for GMJMCMC (Genetically Modified MJMCMC).
#'
#' This is also done when running the algorithm, but this function allows for it to be done manually.
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

#' Function for Calculating Marginal Inclusion Probabilities of Features Given a List of Models
#' @param models The list of models to use.
#'
#' @return A numeric vector of marginal model probabilities based on relative frequencies of model visits in MCMC.
#'
#' @examples
#' result <- gmjmcmc(x = matrix(rnorm(600), 100),
#' y = matrix(rnorm(100), 100), 
#' P = 2, 
#' transforms = c("p0", "exp_dbl"))
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

#' Function for Calculating Feature Importance Through Renormalized Model Estimates
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
  precalc <- matrix(NA, nrow(data$x), length(features))
  for (f in seq_along(features)) {
    feature_string <- print.feature(features[[f]], dataset = TRUE, fixed = data$fixed)
    precalc[, f] <- eval(parse(text = feature_string))
  }
  # Replace any -Inf and Inf values caused by under- or overflow
  precalc <- replace.infinite.data.frame(precalc)
  data$x <- cbind(data$x[, seq_len(data$fixed)], precalc)
  return(data)
}

# TODO: Compare to previous mliks here instead, also add a flag to do that in full likelihood estimation scenarios.
# Function to call the model function
loglik.pre <- function (loglik.pi, model, complex, data, params = NULL, visited.models, sub) {
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
  model.res <- loglik.pi(data$y, data$x, c(rep(TRUE, data$fixed), model), complex, params)
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
check.data <- function (x, y, fixed, verbose) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
    #if (verbose) cat("Data (x) coerced to matrix type.\n")
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
    #if (verbose) cat("Data (y) coerced to matrix type.\n")
  }
  if (nrow(x) != nrow(y)) {
    stop("x and y must have the same number of rows")
  }

  # Ensure that the first F0.size * 2 lines do not contain zero variance variables
  if ((ncol(x) - fixed) * 2 < nrow(x)) {
    vars <- diag(var(x[seq_len((ncol(x) - fixed) * 2), ]))
    for (i in which(vars == 0)[-seq_len(fixed)]) {
      j <- which(x[, i] != x[1, i])[1]
      if (is.na(j)) {
        stop(paste0("column with index ", i, " is constant and only the intercept may be constant, please remove it and try again."))
      }
      x <- rbind(x[j, , drop = FALSE], x[-j, , drop = FALSE])
      y <- rbind(y[j, , drop = FALSE], y[-j, , drop = FALSE])
    }
  }

  return(list(x = x, y = y, fixed = fixed))
}

# Function to extract column names if they are well formed
get.labels <- function (data, verbose) {
  labels <- colnames(data$x)
  if (is.null(labels)) return(FALSE)
  if (sum(is.na(labels)) != 0) {
    if (verbose) cat("NA labels present, using x#\n")
    return(FALSE)
  }
  if (data$fixed > 0) labels <- labels[-seq_len(data$fixed)]
  return(labels)
}
