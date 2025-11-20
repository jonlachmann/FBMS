# Title     : Results functions
# Objective : Functions to handle the results produced
# Created by: jonlachmann
# Created on: 2021-05-06

#' Merge a List of Multiple Results from Many Runs
#'
#' This function will weight the features based on the best marginal posterior in that population
#' and merge the results together, simplifying by merging equivalent features (having high correlation).
#'
#' @param results A list containing multiple results from GMJMCMC (Genetically Modified MJMCMC).
#' @param populations Which populations should be merged from the results, can be "all", "last" (default) or "best".
#' @param complex.measure The complex measure to use when finding the simplest equivalent feature,
#' 1=total width, 2=operation count and 3=depth.
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.0000001
#' @param data Data to use when comparing features, default is NULL meaning that mock data will be generated,
#' if data is supplied it should be of the same form as is required by gmjmcmc, i.e. with both x, y and an intercept.
#'
#' @return An object of class "gmjmcmc_merged" containing the following elements:
#' \item{features}{The features where equivalent features are represented in their simplest form.}
#' \item{marg.probs}{Importance of features.}
#' \item{counts}{Counts of how many versions that were present of each feature.}
#' \item{results}{Results as they were passed to the function.}
#' \item{pop.best}{The population in the results which contained the model with the highest log marginal posterior.}
#' \item{thread.best}{The thread in the results which contained the model with the highest log marginal posterior.}
#' \item{crit.best}{The highest log marginal posterior for any model in the results.}
#' \item{reported}{The highest log marginal likelihood for the reported populations as defined in the populations argument.}
#' \item{rep.pop}{The index of the population which contains reported.}
#' \item{best.log.posteriors}{A matrix where the first column contains the population indices and the second column contains the model with the highest log marginal posterior within that population.}
#' \item{rep.thread}{The index of the thread which contains reported.}
#'
#' @examples
#' result <-  fbms(semimajoraxis ~ ., data = exoplanet,
#'  method = "gmjmcmc.parallel", transforms = c("sigmoid"), 
#'  runs = 2, cores = 1)
#' 
#' summary(result)
#' 
#' plot(result)
#' 
#' merge_results(result$results.raw)
#'
#' @export merge_results
merge_results <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL) {
  # Default values
  if (is.null(populations))
    populations <-"best"
  if (is.null(complex.measure))
    complex.measure <- 2
  if (is.null(tol))
    tol <- 0.0000001
  
  # Check and filter results that did not run successfully
  results <- filter.results(results)
  raw.results <- results
  res.count <- length(results)
  
  # Select populations to use
  res.lengths <- vector("list")
  for (i in 1:res.count) res.lengths[[i]] <- length(results[[i]]$populations)
  if (populations == "last") pops.use <- res.lengths
  else if (populations == "all") pops.use <- lapply(res.lengths, function(x) 1:x)
  else if (populations == "best") pops.use <- lapply(1:res.count, function(x) which.max(unlist(results[[x]]$best.marg)))
  
  # Get the population weigths to be able to weight the features
  pw <- population.weigths(results, pops.use)
  pop.weights <- pw$weights
  
  bests <- matrix(data = 0, ncol = length(results), nrow = length(results[[1]]$populations))
  crit.best <- -Inf
  pop.best <- 1
  thread.best <- 1
  for (i in seq_along(results)) {
    for (pop in 1:(length(results[[i]]$populations))) {
      bests[pop, i] <- results[[i]]$best.margs[[pop]]
      if (results[[i]]$best.margs[[pop]] > crit.best) {
        crit.best <- results[[i]]$best.margs[[pop]]
        pop.best <- pop
        thread.best <- i
      }
    }
  }
  
  # Collect all features and their renormalized weighted values
  features <- vector("list")
  renorms <- vector("list")
  weight_idx <- 1
  for (i in 1:res.count) {
    results[[i]]$pop.weights <- rep(NA, length(results[[i]]$populations))
    results[[i]]$model.probs <- list()
    for (pop in pops.use[[i]]) {
      features <- append(features, results[[i]]$populations[[pop]])
      renorms <- append(renorms, pop.weights[weight_idx] * results[[i]]$marg.probs[[pop]])
      results[[i]]$pop.weights[pop] <- pop.weights[weight_idx]
      weight_idx <- weight_idx + 1
      
      model.probs <- marginal.probs.renorm(results[[i]]$models[[pop]], "models")
      results[[i]]$model.probs[[pop]] <- model.probs$probs
      results[[i]]$models[[pop]] <- results[[i]]$models[[pop]][model.probs$idx]
    }
    accept.tot <- results[[i]]$accept.tot
    best <- results[[i]]$best
    for (item in names(results[[i]])) {
      if (!(item %in% (c("accept.tot", "best", "transforms", "fixed", "intercept", "ncov")))) {
        results[[i]][[item]] <- results[[i]][[item]][pops.use[[i]]]
      }
    }
    results[[i]]$accept.tot <- accept.tot
    results[[i]]$best <- best
  }
  renorms <- unlist(renorms)
  na.feats <- which(is.na(renorms))
  if (length(na.feats) != 0) {
    warning("Underflow occurred,", length(na.feats), "features removed.\n")
    renorms <- renorms[-na.feats]
    features <- features[-na.feats]
  }
  feat.count <- length(features)
  
  # Get complexity for all features
  complex <- complex.features(features)
  
  ## Detect equivalent features
  # Generate mock data to compare features with
  uk <- 1 
  good.mock <- FALSE
  while(!good.mock & uk < 5)
  {
    uk <- uk + 1
    if (is.null(data)) mock.data <- list(x = matrix(runif((results[[1]]$ncov)^2, -100, 100), ncol = results[[1]]$ncov))
    else 
      if(is.null(data$x)) mock.data <- list(x = data) else mock.data = data
    mock.data$fixed = results[[1]]$fixed
    if (results[[1]]$intercept) mock.data$x <- cbind(1, mock.data$x)
    
    mock.data.precalc <- precalc.features(mock.data, features)$x[ , seq_len(feat.count) + results[[1]]$fixed, drop = FALSE]
    
    if(min(sapply(1:dim(mock.data.precalc)[2], function(x)sd(mock.data.precalc[,x])))>0)
    {
      good.mock <- TRUE
      break
    }
  }
 
  cors <- suppressWarnings(cor(mock.data.precalc))
  
  # A map to link equivalent features together,
  # row 1-3 are the simplest equivalent features based on three different complexity measures
  # row 4 is the total weighted density of those features
  feats.map <- matrix(1:feat.count, 4, feat.count, byrow = TRUE)
  for (i in seq_len(nrow(cors))) {
    equiv.feats <- which(cors[i, ] >= (1 - tol))
    # Compare equivalent features complexity to find most simple
    equiv.complex <- list(width = complex$width[equiv.feats], oc = complex$oc[equiv.feats], depth = complex$depth[equiv.feats])
    equiv.simplest <- lapply(equiv.complex, which.min)
    feats.map[1:3, equiv.feats] <- c(equiv.feats[equiv.simplest$width], equiv.feats[equiv.simplest$oc], equiv.feats[equiv.simplest$depth])
    feats.map[4, equiv.feats] <- sum(renorms[equiv.feats])
  }
  # Select the simplest features based on the specified complexity measure and sort them
  feats.simplest.ids <- unique(feats.map[complex.measure, ])
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4, feats.simplest.ids])]
  counts <- sapply(feats.simplest.ids, function(x) sum(feats.map[complex.measure,] == x))
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4, feats.simplest.ids, drop = FALSE]
  merged <- list(
    features = feats.simplest,
    marg.probs = importance,
    counts = counts,
    results = results,
    results.raw = raw.results,
    pop.best = pop.best,
    thread.best = thread.best,
    crit.best = crit.best,
    reported = pw$best,
    rep.pop = pw$pop.best,
    best.log.posteriors = bests,
    rep.thread = pw$thread.best,
    transforms = results[[1]]$transforms,
    fixed = results[[1]]$fixed,
    intercept = results[[1]]$intercept
  )
  attr(merged, "class") <- "gmjmcmc_merged"
  return(merged)
}

filter.results <- function (results) {
  res.count <- length(results)
  res.converged <- sum(sapply(results, function(x) length(x) > 1))
  
  if (res.converged == 0) {
    stop(paste0("All chains resulted in an error!", results[[1]],"\n Please debug and restart"))
  }
  if (res.converged < res.count) {
    warning(paste0("Warning! Some chains resulted in an error: ", results[[which(!sapply(results,function(x)length(x)>1))[1]]],  "'\n Only ",res.converged, " chains finished! \n Only finished chains will be used further!"))
    results <- lapply(results, function (x) {
      if (length(x) > 1) return(x)
      else return(NULL)
    })
    results <- results[sapply(results, function (x) !is.null(x))]
  }
  return(results)
}

# Function for calculating the weights of different populations based on best mlik
population.weigths <- function (results, pops.use) {
  max.crits <- vector("list")
  max.crit <- -Inf
  pop.best <- 1
  thread.best <- 1
  for (i in seq_along(results)) {
    for (pop in pops.use[[i]]) {
      max.crits <- append(max.crits, results[[i]]$best.margs[[pop]])
      if (results[[i]]$best.margs[[pop]] > max.crit) {
        max.crit <- results[[i]]$best.margs[[pop]]
        pop.best <- pop
        thread.best <- i
      }
    }
  }
  max.crits <- unlist(max.crits)
  
  return(list(weights = exp(max.crits-max.crit) / sum(exp(max.crits-max.crit)), best = max.crit, thread.best = thread.best, pop.best = pop.best))
}

#' Function to Generate a Function String for a Model Consisting of Features
#'
#' @param model A logical vector indicating which features to include
#' @param features The population of features
#' @param link The link function to use, as a string
#' @param round Rounding error for the features in the printed format
#'
#' @return A character representation of a model
#'
#' @examples
#' result <- gmjmcmc(y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' P = 2, transforms =  c("p0", "exp_dbl"))
#' summary(result)
#' plot(result)
#' model.string(c(TRUE, FALSE, TRUE, FALSE, TRUE), result$populations[[1]])
#' model.string(result$models[[1]][1][[1]]$model, result$populations[[1]])
#'
#' @export model.string
model.string <- function (model, features, link = "I", round = 2) {
  modelstring <- paste0(sapply(features[model], print.feature, alphas = TRUE, round = round), collapse = "+")
  modelfun <- set_alphas(modelstring)
  modelfun$formula <- paste0(link, "(", modelfun$formula, ")")
  return(modelfun)
}

#' Retrieve the Median Probability Model (MPM)
#'
#' This function extracts the Median Probability Model (MPM) from a fitted model object.
#' The MPM includes features with marginal posterior inclusion probabilities greater than 0.5.
#' It constructs the corresponding model matrix and computes the model fit using the specified likelihood.
#'
#' @param result A fitted model object (e.g., from \code{mjmcmc}, \code{gmjmcmc}, or related classes) containing the summary statistics and marginal probabilities.
#' @param y A numeric vector of response values. For \code{family = "binomial"}, it should contain binary (0/1) responses.
#' @param x A \code{data.frame} of predictor variables. Columns must correspond to features considered during model fitting.
#' @param labels If specified, custom labels of covariates can be used. Default is \code{FALSE}.
#' @param family Character string specifying the model family. Supported options are:
#'   \itemize{
#'     \item \code{"gaussian"} (default) - for continuous outcomes.
#'     \item \code{"binomial"} - for binary outcomes.
#'     \item \code{"custom"} - for user-defined likelihood functions.
#'   }
#' If an unsupported family is provided, a warning is issued and the Gaussian likelihood is used by default.
#' @param loglik.pi A function that computes the log-likelihood. Defaults to \code{gaussian.loglik} unless \code{family = "binomial"}, in which case \code{logistic.loglik} is used. for custom family the user must specify the same likelihood that was used in the inference.
#' @param params Parameters of `loglik.pi`, if not specified NULL will be used by default
#'
#' @return A \code{bgnlm_model} object containing:
#' \describe{
#'   \item{\code{prob}}{The log marginal likelihood of the MPM.}
#'   \item{\code{model}}{A logical vector indicating included features.}
#'   \item{\code{crit}}{Criterion label set to \code{"MPM"}.}
#'   \item{\code{coefs}}{A named numeric vector of model coefficients, including the intercept.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(42)
#' x <- data.frame(
#'   PlanetaryMassJpt = rnorm(100),
#'   RadiusJpt = rnorm(100),
#'   PeriodDays = rnorm(100)
#' )
#' y <- 1 + 0.5 * x$PlanetaryMassJpt - 0.3 * x$RadiusJpt + rnorm(100)
#'
#' # Assume 'result' is a fitted object from gmjmcmc or mjmcmc
#' result <- mjmcmc(cbind(y,x))  
#'
#' # Get the MPM
#' mpm_model <- get.mpm.model(result, y, x, family = "gaussian")
#'
#' # Access coefficients
#' mpm_model$coefs
#' }
#'
#' @export
get.mpm.model <- function(result, y, x, labels = F, family = "gaussian", loglik.pi = gaussian.loglik, params = NULL) {
  
  transforms.bak <- set.transforms(result$transforms)
  
  if (!family %in% c("custom","binomial","gaussian"))
    warning("Unknown family specified. The default gaussian.loglik will be used.")
  
  if(!labels & length(result$labels)>0)
    labels <- result$labels
  
  if (!is.null(attr(result, which = "imputed")))
    x <- impute_x(result,x)
  
  if (family == "binomial")
    loglik.pi <- logistic.loglik
  
  if (is(result, "mjmcmc_parallel")) {
    models <- unlist(lapply(result$chains, function (x) x$models), recursive = FALSE)
    marg.probs <- marginal.probs.renorm(models)$probs
    features <- result$chains[[1]]$populations
  } else if (is(result, "gmjmcmc")) {
    best_pop <- which.max(unlist(result$best.margs))
    marg.probs <- result$marg.probs[[best_pop]]
    features <- result$populations[[best_pop]]
  } else if (is(result, "gmjmcmc_merged")) {
    marg.probs <- result$marg.probs
    features <- result$features
  }else
  {
    marg.probs <- result$marg.probs
    features <- result$populations
  }
  features <- features[marg.probs > 0.5]
  
  if (result$intercept) {
    x <- cbind(1, x)
  }
  precalc <- precalc.features(list(x = x, y = y, fixed = result$fixed), features)
  
  coefs <- loglik.pi(y = y, x = precalc$x, model = rep(TRUE, length(features) + result$fixed), complex = list(oc = 0), mlpost_params = params)$coefs
  
  coefs[is.na(coefs)] <- 0
  
  names(coefs) <- c(names(coefs)[seq_len(result$fixed)], sapply(features, print.feature,labels = labels))
  
  
  model <- structure(list(
    coefs = coefs,
    features = features,
    fixed = result$fixed,
    intercept = result$intercept,
    needs.precalc = FALSE
  ), class = "bgnlm_model")
  
  set.transforms(transforms.bak)
  
  attr(model, which = "imputed") <- attr(result, which = "imputed")
  
  return(model)
}

#' Extract the Best Model from MJMCMC or GMJMCMC Results
#'
#' Retrieves the best model from the results of MJMCMC, MJMCMC parallel, GMJMCMC, or GMJMCMC merged runs 
#' based on the maximum criterion value (\code{crit}). The returned list includes the model probability, 
#' selected features, criterion value, intercept parameter, and named coefficients.
#'
#' @param result An object of class \code{mjmcmc}, \code{mjmcmc_parallel}, \code{gmjmcmc}, or \code{gmjmcmc_merged}, 
#' containing the results from the corresponding model search algorithms.
#' @param labels Logical; if \code{TRUE}, uses labeled feature names when naming the model coefficients. 
#' Alternatively, a character vector of feature names. Default is \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list containing the details of the best model:
#' \describe{
#'   \item{\code{prob}}{A numeric value representing the model's probability.}
#'   \item{\code{model}}{A logical vector indicating which features are included in the best model.}
#'   \item{\code{crit}}{The criterion value used for model selection (e.g., marginal likelihood or posterior probability).}
#'   \item{\code{alpha}}{The intercept parameter of the best model.}
#'   \item{\code{coefs}}{A named numeric vector of model coefficients, including the intercept and selected features.}
#' }
#'
#' @details 
#' The function identifies the best model by selecting the one with the highest \code{crit} value. 
#' Selection logic depends on the class of the \code{result} object:
#' \describe{
#'   \item{\code{mjmcmc}}{Selects the top model from a single MJMCMC run.}
#'   \item{\code{mjmcmc_parallel}}{Identifies the best chain, then selects the best model from that chain.}
#'   \item{\code{gmjmcmc}}{Selects the best population and model within that population.}
#'   \item{\code{gmjmcmc_merged}}{Finds the best chain and population before extracting the top model.}
#' }
#'
#' @examples
#' data(exoplanet)
#' result <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc")
#' get.best.model(result)
#'
#' @export
get.best.model <- function(result, labels = FALSE, ...) {
  if (!inherits(result, c("mjmcmc", "mjmcmc_parallel", "gmjmcmc", "gmjmcmc_merged"))) {
    stop("result must be of class 'mjmcmc', 'mjmcmc_parallel', 'gmjmcmc', or 'gmjmcmc_merged'")
  }
  if (!is.logical(labels) && !is.character(labels)) {
    stop("labels must be a logical value or a character vector")
  }
  UseMethod("get.best.model")
}

#' @method get.best.model mjmcmc
#' @export
get.best.model.mjmcmc <- function(result, labels = FALSE, ...) {
  # Input validation
  stopifnot(
    "result must be of class 'mjmcmc'" = inherits(result, "mjmcmc"),
    "result must contain models" = !is.null(result$models) && length(result$models) > 0,
    "result must contain populations" = !is.null(result$populations) && length(result$populations) > 0
  )
  if (is.character(labels)) {
    stopifnot("labels length must match number of features" = length(labels) == length(result$populations))
  }
  
  # Use result$labels if available and labels is FALSE
  if (length(labels) == 1 && is.logical(labels) && labels[1] == FALSE && length(result$labels) > 0) {
    labels <- result$labels
  }
  
  # Get best model
  best.mod.id <- which.max(sapply(result$models, function(x) x$crit))
  ret <- result$models[[best.mod.id]]
  
  # Assign coefficient names
  coefnames <- sapply(result$populations, print.feature, labels = labels)[ret$model]
  if (result$intercept) coefnames <- c("Intercept", coefnames)
  names(ret$coefs) <- coefnames
  
  # Set attributes and class
  ret$needs.precalc <- FALSE
  class(ret) <- "bgnlm_model"
  attr(ret, which = "imputed") <- attr(result, which = "imputed")
  
  return(ret)
}

#' @method get.best.model mjmcmc_parallel
#' @export
get.best.model.mjmcmc_parallel <- function(result, labels = FALSE, ...) {
  # Input validation
  stopifnot(
    "result must be of class 'mjmcmc_parallel'" = inherits(result, "mjmcmc_parallel"),
    "result must contain chains" = !is.null(result$chains) && length(result$chains) > 0
  )
  if (is.character(labels)) {
    stopifnot("labels length must match number of features" = length(labels) == length(result$chains[[1]]$populations))
  }
  
  # Use result$chains[[1]]$labels if available and labels is FALSE
  if (length(labels) == 1 && is.logical(labels) && labels[1] == FALSE && length(result$chains[[1]]$labels) > 0) {
    labels <- result$chains[[1]]$labels
  }
  
  # Get best chain and model
  best.chain <- which.max(sapply(result$chains, function(x) x$best.crit))
  mod <- get.best.model.mjmcmc(result$chains[[best.chain]], labels = labels)
  attr(mod, which = "imputed") <- attr(result, which = "imputed")
  
  return(mod)
}

#' @method get.best.model gmjmcmc
#' @export
get.best.model.gmjmcmc <- function(result, labels = FALSE, ...) {
  # Input validation
  stopifnot(
    "result must be of class 'gmjmcmc'" = inherits(result, "gmjmcmc"),
    "result must contain models" = !is.null(result$models) && length(result$models) > 0,
    "result must contain populations" = !is.null(result$populations) && length(result$populations) > 0,
    "result must contain best.margs" = !is.null(result$best.margs) && length(result$best.margs) > 0
  )
  if (is.character(labels)) {
    stopifnot("labels length must match number of features" = length(labels) == length(result$populations[[1]]))
  }
  
  # Preserve transforms
  transforms.bak <- set.transforms(result$transforms)
  
  # Use result$labels if available and labels is FALSE
  if (length(labels) == 1 && is.logical(labels) && labels[1] == FALSE && length(result$labels) > 0) {
    labels <- result$labels
  }
  
  # Get best model
  best.pop.id <- which.max(sapply(result$best.margs, function(x) x))
  best.mod.id <- which.max(sapply(result$models[[best.pop.id]], function(x) x$crit))
  ret <- result$models[[best.pop.id]][[best.mod.id]]
  
  # Assign additional attributes
  ret$intercept <- result$intercept
  ret$fixed <- result$fixed
  coefnames <- sapply(result$populations[[best.pop.id]], print.feature, labels = labels)[ret$model]
  if (result$intercept) coefnames <- c("Intercept", coefnames)
  names(ret$coefs) <- coefnames
  
  # Set attributes and class
  ret$needs.precalc <- FALSE
  class(ret) <- "bgnlm_model"
  attr(ret, which = "imputed") <- attr(result, which = "imputed")
  
  # Restore transforms
  set.transforms(transforms.bak)
  
  return(ret)
}

#' @method get.best.model gmjmcmc_merged
#' @export
get.best.model.gmjmcmc_merged <- function(result, labels = FALSE, ...) {
  # Input validation
  stopifnot(
    "result must be of class 'gmjmcmc_merged'" = inherits(result, "gmjmcmc_merged"),
    "result must contain results" = !is.null(result$results) && length(result$results) > 0,
    "result must contain results.raw" = !is.null(result$results.raw) && length(result$results.raw) > 0
  )
  if (is.character(labels)) {
    stopifnot("labels length must match number of features" = length(labels) == length(result$results.raw[[1]]$populations[[1]]))
  }
  
  # Use result$results.raw[[1]]$labels if available and labels is FALSE
  if (length(labels) == 1 && is.logical(labels) && labels[1] == FALSE && length(result$results.raw[[1]]$labels) > 0) {
    labels <- result$results.raw[[1]]$labels
  }
  
  # Get best chain and model
  best.chain <- which.max(sapply(result$results, function(x) x$best))
  mod <- get.best.model.gmjmcmc(result$results.raw[[best.chain]], labels = labels)
  attr(mod, which = "imputed") <- attr(result, which = "imputed")
  
  return(mod)
}


#' Function to Get a Character Representation of a List of Features
#'
#' @param x A list of feature objects
#' @param round Rounding precision for parameters of the features
#'
#' @return A matrix of character representations of the features of a model.
#'
#' @examples
#' result <- gmjmcmc(y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' P = 2, 
#' transforms = c("p0", "exp_dbl"))
#' string.population(result$populations[[1]])
#'
#' @export
string.population <- function(x, round = 2) {
  cbind(sapply(x, print.feature, round = round))
}

#' Function to Get a Character Representation of a List of Models
#'
#' @param features A list of feature objects on which the models are build
#' @param models A list of model objects
#' @param round Rounding precision for parameters of the features
#' @param link The link function to use, as a string
#'
#' @return A matrix of character representations of a list of models.
#'
#' @examples
#' result <- gmjmcmc(y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' P = 2, 
#' transforms = c("p0", "exp_dbl"))
#' string.population.models(result$populations[[2]], result$models[[2]])
#'
#' @export
string.population.models <- function(features, models, round = 2, link = "I") {
  cbind(sapply(seq_along(models), FUN = function(x) model.string(features = features, model = (models[[x]]$model), round = round, link = "I")))
}

#' Function to Plot GMJMCMC Results and Merged Results from merge.results
#'
#' @param x The results to use
#' @param count The number of features to plot, defaults to all
#' @param pop The population to plot, defaults to last
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.0000001
#' @param data Data to merge on, important if pre-filtering was used
#' @param ... Not used.
#'
#' @return No return value, just creates a plot
#'
#' @examples
#' result <- gmjmcmc(y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' P = 2, 
#' transforms = c("p0", "exp_dbl"))
#' plot(result)
#' 
#'
#' @export 
plot.gmjmcmc <- function (x, count = "all", pop = "best", tol = 0.0000001, data = NULL, ...) {
  transforms.bak <- set.transforms(x$transforms)
  if (pop != "last") {
    results <- list()
    results[[1]] <- x
    x <- merge_results(results, pop, 2, 0.0000001, data = data)
    return(marg.prob.plot(sapply(x$features, print), x$marg.probs, count = count))
  }
  
  if (pop == "last") pop <- length(x$populations)
  if (is.null(x$populations)) {
    pops <- x$features
    marg.probs <- x$marg.probs
  } else {
    pops <- x$populations[[pop]]
    marg.probs <- x$marg.probs[[pop]]
  }
  plot.mjmcmc(list(populations = pops, marg.probs = marg.probs), count)
  set.transforms(transforms.bak)
  return("done")
}


#' Plot BGNLM Model
#'
#' Plots the coefficients of a BGNLM model.
#'
#' @param x Object of class "bgnlm_model".
#' @param ... Additional arguments passed to barplot.
#' @return The input object (invisibly).
#' @method plot bgnlm_model
#' @importFrom stats coef
#' @export
#' @examples
#' data(exoplanet)
#' model <- get.best.model(fbms(semimajoraxis ~ ., data = exoplanet, family = "gaussian"))
#' plot(model)
plot.bgnlm_model <- function(x, ...) {
  stopifnot(inherits(x, "bgnlm_model"))
  coefs <- coef(x)
  barplot(coefs, main = "BGNLM Coefficients", names.arg = names(coefs), horiz  = TRUE, ...)
}

#' Function to Plot GMJMCMC Results and Merged Results from merge.results
#'
#' @param x The results to use
#' @param count The number of features to plot, defaults to all
#' @param ... Not used.
#'
#' @return No return value, just creates a plot
#'
#' @examples
#' result <- mjmcmc(
#' y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100),
#' loglik.pi = gaussian.loglik)
#' plot(result)
#'
#' @export 
plot.mjmcmc <- function (x, count = "all", ...) {
  transforms.bak <- set.transforms(x$transforms)
  ## Get features as strings for printing and marginal posteriors
  # If this is a merged results the structure is one way
  if (is.null(x$populations)) {
    feats.strings <- sapply(x$features, print)
    feats.strings <- paste0(feats.strings, ", ", x$count)
    marg.probs <- x$marg.probs
  } # If this is a result that is not merged, it is another way
  else {
    feats.strings <- sapply(x$populations, print)
    marg.probs <- x$marg.probs
  }
  
  marg.prob.plot(feats.strings, marg.probs, count)
  set.transforms(transforms.bak)
  return("done")
}

marg.prob.plot <- function (feats.strings, marg.probs, count = "all", ...) {
  # Plot the distribution
  feats.strings <- feats.strings[order(marg.probs)]
  marg.probs <- sort(marg.probs)
  tot <- length(marg.probs)
  if (count == "all") count <- tot
  y <- barplot(marg.probs[(tot - count + 1):tot], horiz = TRUE, xlab = "Marginal probability", ylab = "Feature")
  text((max(marg.probs[(tot - count + 1):tot]) / 2), y, feats.strings[(tot - count + 1):tot])
}

#' Plot an mjmcmc_parallel Run
#' @inheritParams plot.mjmcmc
#' @return No return value, just creates a plot
#' 
#' @examples
#' result <- mjmcmc.parallel(runs = 1, 
#' cores = 1, 
#' y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' loglik.pi = gaussian.loglik)
#' plot(result)
#' 
#' @export 
plot.mjmcmc_parallel <- function (x, count = "all", ...) {
  merged <- merge_mjmcmc_parallel(x)
  marg.prob.plot(merged$features, merged$marg.probs, count)
}

merge_mjmcmc_parallel <- function (x) {
  run.weights <- run.weigths(x)
  marg.probs <- x$chains[[1]]$marg.probs * run.weights[1]
  for (i in seq_along(x[-c(1, (-1:0 + length(x)))])) {
    marg.probs <- marg.probs + x$chains[[i]]$marg.probs * run.weights[i]
  }
  return(structure(
    list(
      features = sapply(x$chains[[1]]$populations, print),
      marg.probs = marg.probs,
      results = x
    ),
    class = "mjmcmc_merged"
  ))
}


run.weigths <- function (results) {
  best.crits <- sapply(results$chains, function (x) x$best.crit)
  max.crit <- max(best.crits)
  return(exp(best.crits - max.crit) / sum(exp(best.crits - max.crit)))
}

#' Plot a gmjmcmc_merged Run
#' @inheritParams plot.gmjmcmc
#' @return No return value, just creates a plot
#' 
#' @examples
#' result <- gmjmcmc.parallel(
#'  runs = 1,
#'  cores = 1,
#'  y = matrix(rnorm(100), 100),
#'  x = matrix(rnorm(600), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' plot(result)
#' 
#' @export 
plot.gmjmcmc_merged <- function (x, count = "all", pop = NULL,tol =  0.0000001, data = NULL, ...) {
  transforms.bak <- set.transforms(x$transforms)
  if (!is.null(pop)) {
    x <- merge_results(x$results.raw, pop, 2, 0.0000001, data = data)
  }
  
  marg.prob.plot(sapply(x$features[x$marg.probs > tol], print), x$marg.probs[x$marg.probs > tol], count = count)
  set.transforms(transforms.bak)
  return("done")
}


#' Compute Effects for Specified Covariates Using a Fitted Model
#'
#' This function computes model averaged effects for specified covariates using a fitted model object.
#' The effects are expected change in the BMA linear predictor having an increase of the corresponding covariate by one unit, while other covariates are fixed to 0.
#' Users can provide custom labels and specify quantiles for the computation of effects.
#'
#' @param object A fitted model object, typically the result of a regression or predictive modeling.
#' @param labels A vector of labels for which effects are to be computed.
#' @param quantiles A numeric vector specifying the quantiles to be calculated. Default is c(0.025, 0.5, 0.975).
#'
#' @return A matrix of treatment effects for the specified labels, with rows corresponding to labels and columns to quantiles.
#'
#' @examples
#'
#' data <- data.frame(matrix(rnorm(600), 100))
#' result <- mjmcmc.parallel(runs = 2, 
#' cores = 1, 
#' y = matrix(rnorm(100), 100),
#' x = data, 
#' loglik.pi = gaussian.loglik)
#' compute_effects(result,labels = names(data))
#'
#' @seealso \code{\link{predict}}
#' @export
compute_effects <- function(object, labels, quantiles = c(0.025, 0.5, 0.975)) {
  effects <- rbind(0, diag(length(labels)))
  preds.eff <- predict(object = object, x = as.matrix(effects), quantiles = quantiles)
  if (length(preds.eff$aggr) > 0)
    preds.eff <- t(preds.eff$aggr$quantiles)
  else
    preds.eff <- t(preds.eff$quantiles)
  preds.eff[2:(length(labels) + 1), ] <- preds.eff[2:(length(labels) + 1), ] - preds.eff[1, ]
  
  summ <- data.frame(cbind(c("intercept", labels), round(preds.eff, 4)))
  names(summ) <- c("Covariate", paste0("quant_", quantiles))
  return(summ)
}
