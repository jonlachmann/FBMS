#' Function to Print a Quick Summary of the Results
#'
#' @param object The results to use
#' @param pop The population to print for, defaults to last
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param labels Should the covariates be named, or just referred to as their place in the data.frame.
#' @param effects Quantiles for posterior modes of the effects across models to be reported, if either effects are NULL or if labels are NULL, no effects are reported.
#' @param data Data to merge on, important if pre-filtering was used
#' @param verbose If the summary should be printed to the console or just returned, defaults to TRUE
#' @param ... Not used.
#'
#' @return A data frame containing the following columns:
#' \item{feats.strings}{Character representation of the features ordered by marginal probabilities.}
#' \item{marg.probs}{Marginal probabilities corresponding to the ordered feature strings.}
#'
#' @examples
#' result <- gmjmcmc(y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' P = 2, 
#' transforms =  c("p0", "exp_dbl"))
#' summary(result, pop = "best")
#'
#' @export
summary.gmjmcmc <- function (object, pop = "best", tol = 0.0001, labels = FALSE, effects = NULL, data = NULL, verbose = TRUE, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if (length(labels) == 1 && labels[1] == FALSE && length(object$labels) > 0) {
    labels = object$labels
  }
  if (pop == "all") {
    results <- list()
    results[[1]] <- object
    merged <- merge_results(results, pop, 2, 0.0000001, data = data)

    best <- max(sapply(merged$results, function (y) y$best))
    feats.strings <- sapply(merged$features, FUN = function(x) print.feature(x = x, labels = labels, round = 2))

    if (!is.null(effects) & !is.null(labels)) {
      effects <- compute_effects(merged,labels = labels, quantiles = effects)
    }

    return(summary_internal(
      best = merged$crit.best,
      feats.strings,
      merged$marg.probs,
      effects = effects,
      best.pop = merged$pop.best,
      thread.best = merged$thread.best,
      reported = merged$reported,
      rep.pop = merged$rep.pop,
      rep.thread = merged$rep.thread,
      tol = tol,
      verbose = verbose
    ))
  }

  if (pop == "last") pop <- length(object$models)
  else if (pop == "best") pop <- which.max(unlist(object$best.margs))
  feats.strings <- sapply(object$populations[[pop]], FUN = function(x) print.feature(x = x, labels = labels, round = 2))

  if (!is.null(effects) & !is.null(labels)) {
    effects <- compute_effects(object, labels = labels, quantiles = effects)
  }

  obj <- summary_internal(
    best = object$best,
    marg.probs = object$marg.probs[[pop]],
    effects = effects,
    feats.strings = feats.strings,
    best.pop = which.max(unlist(object$best.margs)),
    reported = object$best.margs[[pop]],
    rep.pop = pop,
    tol = tol,
    verbose = verbose
  )
  set.transforms(transforms.bak)
  return(obj)
}



#' Function to Print a Quick Summary of the Results
#'
#' @param object The results to use
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param labels Should the covariates be named, or just referred to as their place in the data.frame.
#' @param effects Quantiles for posterior modes of the effects across models to be reported, if either effects are NULL or if labels are NULL, no effects are reported.
#' @param pop If null same as in merge.options for running parallel gmjmcmc otherwise results will be re-merged according to pop that can be "all", "last", "best"
#' @param data Data to merge on, important if pre-filtering was used
#' @param verbose If the summary should be printed to the console or just returned, defaults to TRUE
#' @param ... Not used.
#'
#' @return A data frame containing the following columns:
#' \item{feats.strings}{Character representation of the features ordered by marginal probabilities.}
#' \item{marg.probs}{Marginal probabilities corresponding to the ordered feature strings.}
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
#' summary(result)
#'
#' @export
summary.gmjmcmc_merged <- function (object, tol = 0.0001, labels = FALSE, effects = NULL, pop = NULL, data = NULL, verbose = TRUE, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if (length(labels) == 1 && labels[1] == FALSE && length(object$results.raw[[1]]$labels) > 0) {
    labels = object$results.raw[[1]]$labels
  }
  if (!is.null(pop)) {
    object <- merge_results(object$results.raw, populations = pop, complex.measure = 2, tol = 0.0000001, data = data)
  }

  best <- max(sapply(object$results, function (y) y$best))
  feats.strings <- sapply(object$features, FUN = function(x) print.feature(x = x, labels = labels, round = 2))

  if (!is.null(effects) & !is.null(labels)) {
    effects <- compute_effects(object,labels = labels, quantiles = effects)
  }

  obj <- summary_internal(
    best = object$crit.best,
    feats.strings = feats.strings,
    marg.probs = object$marg.probs,
    effects = effects,
    best.pop = object$pop.best,
    thread.best = object$thread.best,
    reported = object$reported,
    rep.pop = object$rep.pop,
    rep.thread = object$rep.thread,
    tol = tol,
    verbose = verbose
  )
  set.transforms(transforms.bak)
  return(obj)
}

#' Function to Print a Quick Summary of the Results
#'
#' @param object The results to use
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param labels Should the covariates be named, or just referred to as their place in the data.frame.
#' @param effects Quantiles for posterior modes of the effects across models to be reported, if either effects are NULL or if labels are NULL, no effects are reported.
#' @param verbose If the summary should be printed to the console or just returned, defaults to TRUE
#' @param ... Not used.
#'
#' @return A data frame containing the following columns:
#' \item{feats.strings}{Character representation of the covariates ordered by marginal probabilities.}
#' \item{marg.probs}{Marginal probabilities corresponding to the ordered feature strings.}
#'
#' @examples
#' result <- mjmcmc(y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' loglik.pi = gaussian.loglik)
#' summary(result)
#'
#' @export
summary.mjmcmc <- function (object, tol = 0.0001, labels = FALSE, effects = NULL, verbose = TRUE, ...) {
  return(summary.mjmcmc_parallel(
    list(chains = list(object), fixed = object$fixed, intercept = object$intercept),
    tol = tol,
    labels = labels,
    effects = effects,
    verbose = verbose
  ))
}

#' Function to Print a Quick Summary of the Results
#'
#' @param object The results to use
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param labels Should the covariates be named, or just referred to as their place in the data.frame.
#' @param effects Quantiles for posterior modes of the effects across models to be reported, if either effects are NULL or if labels are NULL, no effects are reported.
#' @param verbose If the summary should be printed to the console or just returned, defaults to TRUE
#' @param ... Not used.
#'
#' @return A data frame containing the following columns:
#' \item{feats.strings}{Character representation of the covariates ordered by marginal probabilities.}
#' \item{marg.probs}{Marginal probabilities corresponding to the ordered feature strings.}
#'
#' @examples
#' result <- mjmcmc.parallel(runs = 1, 
#' cores = 1,  
#' y = matrix(rnorm(100), 100),
#' x = matrix(rnorm(600), 100), 
#' loglik.pi = gaussian.loglik)
#' summary(result)
#'
#' @export
summary.mjmcmc_parallel <- function (object, tol = 0.0001, labels = FALSE, effects = NULL, verbose = TRUE, ...) {
  # Get features as strings for printing
  if (length(labels) == 1 && labels[1] == FALSE && length(object$chains[[1]]$labels) > 0) {
    labels = object$chains[[1]]$labels
  }
  feats.strings <- sapply(object$chains[[1]]$populations, FUN = function(x) print.feature(x = x, labels = labels, round = 2))
  # Get marginal posterior of features

  models <- unlist(lapply(object$chains, function (x) x$models), recursive = FALSE)
  marg.probs <- marginal.probs.renorm(models)$probs
  best <- max(sapply(object$chains, function (x) x$best))
  if (!is.null(effects) & !is.null(labels)) {
    effects <- compute_effects(object$chains[[1]], labels = labels, quantiles = effects)
  }
  return(summary_internal(best, feats.strings, marg.probs, effects, tol = tol, verbose = verbose))
}

summary_internal <- function (best, feats.strings, marg.probs, effects = NULL, tol = 0.0001, best.pop = NULL, reported = NULL, rep.pop = NULL, rep.thread = NULL, thread.best = NULL, verbose = TRUE) {
  keep <- which(marg.probs[1, ] > tol)

  if (verbose) {
    # Print the best log marginal posterior
    if (length(best.pop) > 0) {
      if (length(thread.best) > 0) {
        cat("\nBest   population:", best.pop, " thread:", thread.best, " log marginal posterior:", best,"\n")
        if(best.pop!=rep.pop)
          cat("Report population:", rep.pop, " thread:", rep.thread, " log marginal posterior:", reported,"\n")
      } else {
        cat("\nBest   population:", best.pop, " log marginal posterior:", best,"\n")
        if(best.pop!=rep.pop)
          cat("Report population:", rep.pop, " log marginal posterior:", reported,"\n")
      }
    } else {
      cat("\nBest log marginal posterior: ", best,"\n")
    }
    cat("\n")
  }

  feats.strings <- feats.strings[keep]
  marg.probs <- marg.probs[1, keep]
  ord.marg <- order(marg.probs, decreasing = TRUE)

  if (!is.null(effects)) {
    return(list(PIP = data.frame(feats.strings = feats.strings[ord.marg], marg.probs = marg.probs[ord.marg]), EFF = effects))
  }

  return(data.frame(feats.strings = feats.strings[ord.marg], marg.probs = marg.probs[ord.marg]))
}
