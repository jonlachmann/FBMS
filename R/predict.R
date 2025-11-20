#' Predict Responses from a BGNLM Model
#'
#' This function generates predictions from a fitted \code{bgnlm_model} object given a new dataset.
#'
#' @param object A fitted \code{bgnlm_model} object obtained from the BGNLM fitting procedure. 
#'              It should contain the estimated coefficients in \code{model$coefs}.
#' @param x A \code{data.frame} containing the new data for which predictions are to be made. 
#'          The variables in \code{x} must match the features used in the model.
#' @param link A link function to apply to the linear predictor. 
#'             By default, it is the identity function \code{function(x)\{x\}}, 
#'             but it can be any function such as \code{plogis} for logistic regression models.
#' @param x_train Training design matrix to be provided when imputations are to be made from them
#' @param ...  Additional arguments to pass to prediction function.
#'
#' @return A numeric vector of predicted values for the given data \code{x}. 
#'         These predictions are calculated as \eqn{\hat{y} = \text{link}(X \beta)}, 
#'         where \eqn{X} is the design matrix and \eqn{\beta} are the model coefficients.
#'
#' @examples
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet)
#' preds <- predict(get.best.model(model), exoplanet[,-1])
#' @export
predict.bgnlm_model <- function(object, x, link = function(x) {x}, x_train = NULL, ... ) {
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x_test = x, x_train = x_train)
  x <- data.frame(x)
  if (object$needs.precalc) {
    if (object$intercept) {
      x <- cbind(1, x)
    }
    
    if(length(object$features)==0)
    {
      warning("MPM has no featres included! All posteriors below 0.5! Baseline only used.")
      x.precalc <-  model.matrix(~1, data = x)
    } else precalc <- precalc.features(list(x = x, y = NULL, fixed = object$fixed), object$features)
    
    if (dim(precalc$x)[2]>length(object$coefs[object$coefs!=0])) {
      precalc$x <- as.matrix(precalc$x[,-1])
    }
    
    yhat <- link(precalc$x %*% object$coefs[object$coefs != 0])
  } else {
    
    if(length(object$coefs)==1)
    {
      warning("MPM has no featres included! All posteriors below 0.5! Baseline only used.")
      x.precalc <-  model.matrix(~1, data = x)
    }
    else{
      x.precalc <- model.matrix(
        as.formula(paste0("~I(", paste0(names(object$coefs)[-1][object$coefs[-1]!=0], collapse = ")+I("), ")")),
        data = x
      )
    }
    
    if (dim(x.precalc)[2]<length(object$coefs[object$coefs!=0])) {
      x.precalc <- cbind(1,x.precalc)
    } else if (dim(x.precalc)[2]>length(object$coefs[object$coefs!=0])) {
      x.precalc <- as.matrix(x.precalc[,-1])
    }
    yhat <- link(x.precalc %*% object$coefs[object$coefs!=0])
  }
  return(yhat)
}


#' Predict Using a GMJMCMC Result Object
#'
#' @inheritParams predict.gmjmcmc_merged
#' @return A list containing aggregated predictions and per model predictions.
#' \item{aggr}{Aggregated predictions with mean and quantiles.}
#' \item{preds}{A list of lists containing individual predictions per model per population in object.}
#' 
#' @examples
#' result <- gmjmcmc(
#'  x = matrix(rnorm(600), 100),
#'  y = matrix(rnorm(100), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' 
#' @export
predict.gmjmcmc <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975),  pop = NULL,tol =  0.0000001, x_train = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  merged <- merge_results(list(object), data = list(x = x, object$fixed), populations = pop, tol = tol)
  set.transforms(transforms.bak)
  return(predict.gmjmcmc_merged(merged, x, link, quantiles))
}

#' New Idea for a More Streamlined Function...
#' Produces slightly different results from the fun above since this is using all lo.models too.
#' @inheritParams predict.gmjmcmc_merged
#' @param pop The population to use.
#' @noRd
predict.gmjmcmc.2 <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = 1, x_train = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  mmodel <- lapply(object[1:8], function (x) x[[pop]])
  
  # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
  x.precalc <- precalc.features(cbind(0, 1, x), mmodel$populations)[, -1]
  set.transforms(transforms.bak)
  return(predict.mjmcmc(mmodel, x.precalc, link, quantiles))
}

#' Predict Using a Merged GMJMCMC Result Object
#'
#' @param object The model to use.
#' @param x The new data to use for the prediction, a matrix where each row is an observation.
#' @param link The link function to use
#' @param quantiles The quantiles to calculate credible intervals for the posterior modes (in model space).
#' @param pop The population to plot, defaults to last
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.0000001
#' @param x_train Training design matrix to be provided when imputations are to be made from them
#' 
#' @param ... Not used.
#' @return A list containing aggregated predictions and per model predictions.
#' \item{aggr}{Aggregated predictions with mean and quantiles.}
#' \item{preds}{A list of lists containing individual predictions per model per population in object.}
#'
#' @examples
#' result <- gmjmcmc.parallel(
#'  runs = 1,
#'  cores = 1,
#'  x = matrix(rnorm(600), 100),
#'  y = matrix(rnorm(100), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' preds <- predict(result, matrix(rnorm(600), 100))
#'
#' @export
predict.gmjmcmc_merged <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = NULL, tol = 0.0000001, x_train = NULL, ...) {
  
  
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  if (object$intercept) {
    x <- cbind(1, x)
  }
  
  
  transforms.bak <- set.transforms(object$transforms)
  if (!is.null(pop))
    object <- merge_results(object$results.raw, pop, 2, tol, data = list(x = x, fixed = object$fixed))
  
  preds <- list()
  for (i in seq_along(object$results)) {
    preds[[i]] <- list()
    for (j in seq_along(object$results[[i]]$populations)) {
      # Select the models and features to predict from at this iteration
      models <- object$results[[i]]$models[[j]]
      features <- object$results[[i]]$populations[[j]]
      model.probs <- object$results[[i]]$model.probs[[j]]
      
      # Precalculate the features for the new data
      x.precalc <- precalc.features(list(x = x, fixed = object$fixed), features)$x
      
      yhat <- matrix(0, nrow = nrow(x), ncol = length(models))
      for (k in seq_along(models)) {
        # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
        if (models[[k]]$crit == -.Machine$double.xmax) next
        yhat[, k] <- link(x.precalc[, c(rep(TRUE, object$fixed), models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
      }
      
      mean.pred <- rowSums(yhat %*% diag(as.numeric(model.probs)))
      pred.quant <- apply(yhat, 1, weighted.quantiles, weights=model.probs, prob=quantiles)
      
      preds[[i]][[j]] <- list(mean=mean.pred, quantiles=pred.quant, weight=object$results[[i]]$pop.weights[j])
    }
  }
  
  aggr <- list()
  aggr$mean <- 0 * preds[[1]][[1]]$mean
  aggr$quantiles <- 0 * preds[[1]][[1]]$quantiles
  for (i in seq_along(preds)) {
    for (j in seq_along(preds[[i]])) {
      aggr$mean <- aggr$mean + preds[[i]][[j]]$mean * object$results[[i]]$pop.weights[j]
      aggr$quantiles <- aggr$quantiles + preds[[i]][[j]]$quantiles * object$results[[i]]$pop.weights[j]
    }
  }
  set.transforms(transforms.bak)
  
  result <- list(aggr = aggr, preds = preds)
  class(result) <- "fbms_predict"
  return(result)
}

#' Predict Using an MJMCMC Result Object
#'
#' @inheritParams predict.gmjmcmc_merged
#' @return A list containing aggregated predictions.
#' \item{mean}{Mean of aggregated predictions.}
#' \item{quantiles}{Quantiles of aggregated predictions.}
#' 
#' @examples
#' result <- mjmcmc(
#' x = matrix(rnorm(600), 100),
#' y = matrix(rnorm(100), 100),
#' loglik.pi = gaussian.loglik)
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' @export
predict.mjmcmc <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), x_train = NULL, ...) {
  # Select the models and features to predict from at this iteration
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  
  if (object$intercept) {
    x <- cbind(1, x)
  }
  
  
  models <- object$models[object$model.probs.idx]
  
  yhat <- matrix(0, nrow = nrow(x), ncol = length(models))
  for (k in seq_along(models)) {
    # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
    if (models[[k]]$crit == -.Machine$double.xmax) next
    yhat[, k] <- link(x[, c(rep(TRUE, object$fixed), models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
  }
  
  mean.pred <- rowSums(yhat %*% diag(as.numeric(object$model.probs)))
  pred.quant <- apply(yhat, 1, weighted.quantiles, weights = object$model.probs, prob = quantiles)
  
  result <- list(aggr = list(mean = mean.pred, quantiles = pred.quant), preds = mean.pred)
  class(result) <- "fbms_predict"
  return(result)
}

#' Predict Using an MJMCMC Result Object from a Parallel Run
#'
#' @inheritParams predict.gmjmcmc_merged
#' @return A list containing aggregated predictions.
#' \item{mean}{Mean of aggregated predictions.}
#' \item{quantiles}{Quantiles of aggregated predictions.}
#' 
#' @examples
#' result <- mjmcmc.parallel(runs = 1, 
#' cores = 1, 
#' x = matrix(rnorm(600), 100),
#' y = matrix(rnorm(100), 100), 
#' loglik.pi = gaussian.loglik)
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' @export
predict.mjmcmc_parallel <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), x_train = NULL, ...) {
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  max.crits <- sapply(object$chains, function (x) x$best.crit)
  max.crit <- max(max.crits)
  result.weights <- exp(max.crits - max.crit) / sum(exp(max.crits - max.crit))
  
  preds <- lapply(object$chains, predict,x, link, quantiles)
  
  aggr <- list()
  aggr$mean <- 0 * preds[[1]]$aggr$mean
  aggr$quantiles <- 0 * preds[[1]]$aggr$quantiles
  for (i in seq_along(preds)) {
    aggr$mean <- aggr$mean + preds[[i]]$aggr$mean * result.weights[i]
    aggr$quantiles <- aggr$quantiles + preds[[i]]$aggr$quantiles * result.weights[i]
  }
  
  result <- list(aggr = aggr, preds = preds)
  class(result) <- "fbms_predict"
  return(result)
  
}

#' Predict Using a GMJMCMC Result Object from a Parallel Run
#'
#' @inheritParams predict.gmjmcmc_merged
#' @param ... Additional arguments to pass to merge_results.
#' @return A list containing aggregated predictions and per model predictions.
#' \item{aggr}{Aggregated predictions with mean and quantiles.}
#' \item{preds}{A list of lists containing individual predictions per model per population in object.}
#' 
#' @examples
#' result <- gmjmcmc.parallel(
#'  runs = 1,
#'  cores = 1,
#'  x = matrix(rnorm(600), 100),
#'  y = matrix(rnorm(100), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' @export
predict.gmjmcmc_parallel <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), x_train = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  merged <- merge_results(object,data = cbind(1, x), ...)
  results <- predict.gmjmcmc_merged(merged, x, link, quantiles)
  set.transforms(transforms.bak)
  return(results)
}

#' Calculate Weighted Quantiles
#'
#' @param values The values to use
#' @param weights The weights of the values
#' @param prob The probabilities of the quantiles to use
#'
#' @return Weighted quantiles
#' @noRd
weighted.quantiles <- function (values, weights, prob = c(0.025, 0.975)) {
  ordered <- order(values)
  P <- cumsum(weights[ordered])
  
  iv <- integer(length(prob))
  for (i in seq_along(iv)) {
    iv[i] <- which.max(P >= prob[i])
  }
  {values[ordered]}[iv]
}



#' Generic for Accessing Aggregated Predictions
#'
#' Dispatches to methods for extracting aggregated predictions from objects.
#'
#' @param object An object.
#' @param ... Additional arguments passed to methods.
#' @return Aggregated predictions (format depends on the object class).
#' @export
aggr <- function(object, ...) UseMethod("aggr")

#' Access Aggregated Predictions
#'
#' Extracts the aggregated predictions (mean and quantiles) from an FBMS prediction object.
#'
#' @param object Object of class "fbms_predict".
#' @param ... Additional arguments (ignored).
#' @return List containing aggregated mean and quantiles.
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet)
#' pred <- predict(model, exoplanet[51:60, -1])
#' aggr(pred)
#' }
aggr.fbms_predict <- function(object, ...) {
  stopifnot(inherits(object, "fbms_predict"))
  object$aggr
}

#' Generic for Accessing Quantile Predictions
#'
#' Dispatches to methods for extracting quantile predictions from objects.
#'
#' @param object An object.
#' @param ... Additional arguments passed to methods.
#' @return Quantile predictions (format depends on the object class).
#' @export
predquantiles <- function(object, ...) UseMethod("predquantiles")

#' Access Quantile Predictions
#'
#' Extracts the quantile predictions from an FBMS prediction object.
#'
#' @param object Object of class "fbms_predict".
#' @param ... Additional arguments (ignored).
#' @return Matrix of quantile predictions, or NULL if not available.
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc")
#' pred <- predict(model, exoplanet[51:60, -1])
#' predquantiles(pred)
#' }
predquantiles.fbms_predict <- function(object, ...) {
  stopifnot(inherits(object, "fbms_predict"))
  object$aggr$quantiles
}


#' Generic for Accessing Mean Predictions
#'
#' Dispatches to methods for extracting quantile predictions from objects.
#'
#' @param object An object.
#' @param ... Additional arguments passed to methods.
#' @return Posterior mean predictions (format depends on the object class).
#' @export
predmean <- function(object, ...) UseMethod("predmean")


#' Access Mean Predictions
#'
#' Extracts the mean predictions from an FBMS prediction object.
#'
#' @param object Object of class "fbms_predict".
#' @param ... Additional arguments (ignored).
#' @return Vector of mean predictions.
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet, method = "mjmcmc")
#' pred <- predict(model, exoplanet[51:60, -1])
#' predmean(pred)
#' }
predmean.fbms_predict <- function(object, ...) {
  stopifnot(inherits(object, "fbms_predict"))
  object$aggr$mean
}

#' Access Fitted Values
#'
#' Extracts the mean predictions from an FBMS prediction object (alias for mean).
#'
#' @param object Object of class "fbms_predict".
#' @param ... Additional arguments (ignored).
#' @return Vector of mean predictions.
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet)
#' pred <- predict(model, exoplanet[51:60, -1])
#' fitted(pred)
#' }
fitted.fbms_predict <- function(object, ...) {
  stopifnot(inherits(object, "fbms_predict"))
  object$aggr$mean
}

#' Print FBMS Prediction Object
#'
#' Displays a summary of an FBMS prediction object, including mean predictions and quantiles.
#'
#' @param x Object of class "fbms_predict".
#' @param ... Additional arguments (ignored).
#' @return Prints a summary and returns NULL.
#' @method print fbms_predict
#' @importFrom utils head
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet)
#' pred <- predict(model, exoplanet[51:60, -1])
#' print(pred)
#' }
print.fbms_predict <- function(x, ...) {
  stopifnot(inherits(x, "fbms_predict"))
  cat("FBMS Prediction Object:\n")
  cat("  Number of predictions:", length(x$aggr$mean), "\n")
  cat("  Mean Predictions (first 6):\n")
  print(head(x$aggr$mean, 6))
  if (!is.null(x$aggr$quantiles)) {
    cat("  Quantiles (first 6 rows):\n")
    print(head(x$aggr$quantiles, 6))
  }
  if (!is.null(x$preds)) {
    cat("  Number of populations:", length(x$preds), "\n")
  }
  invisible(NULL)
}

#' Summary of FBMS Prediction Object
#'
#' Provides a detailed summary of an FBMS prediction object, including prediction ranges.
#'
#' @param object Object of class "fbms_predict".
#' @param ... Additional arguments (ignored).
#' @return Prints a summary and returns NULL.
#' @method summary fbms_predict
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet)
#' pred <- predict(model, exoplanet[51:60, -1])
#' summary(pred)
#' }
summary.fbms_predict <- function(object, ...) {
  stopifnot(inherits(object, "fbms_predict"))
  cat("Summary of FBMS Predictions:\n")
  cat("  Number of predictions:", length(object$aggr$mean), "\n")
  cat("  Mean prediction range:", format(range(object$aggr$mean), digits = 3), "\n")
  if (!is.null(object$aggr$quantiles)) {
    cat("  Quantile ranges:\n")
    print(apply(object$aggr$quantiles, 2, function(x) format(range(x), digits = 3)))
  }
  if (!is.null(object$preds)) {
    cat("  Number of populations:", length(object$preds), "\n")
  }
  invisible(NULL)
}

#' Plot FBMS Prediction Object
#'
#' Plots the mean predictions and quantile intervals from an FBMS prediction object, with quantiles in varying shades of grey.
#'
#' @param x Object of class "fbms_predict".
#' @param ... Additional arguments passed to plot.
#' @return Plots the predictions and returns NULL.
#' @method plot fbms_predict
#' @importFrom grDevices grey.colors
#' @importFrom graphics legend matlines
#' @export
#' @examples
#' \donttest{
#' data(exoplanet)
#' model <- fbms(semimajoraxis ~ ., data = exoplanet)
#' pred <- predict(model, exoplanet[51:60, -1], 
#' quantiles = c(0.025, 0.1, 0.5, 0.9, 0.975))
#' plot(pred)
#' }
plot.fbms_predict <- function(x, ...) {
  stopifnot(inherits(x, "fbms_predict"))
  n <- length(x$aggr$mean)
  
  # Validate quantile dimensions
  if (!is.null(x$aggr$quantiles)) {
    if (!is.matrix(x$aggr$quantiles)) {
      x$aggr$quantiles <- as.matrix(x$aggr$quantiles)
    }
    if (nrow(x$aggr$quantiles) != n) {
      if (ncol(x$aggr$quantiles) == n) {
        x$aggr$quantiles <- t(x$aggr$quantiles)  # Transpose if needed
      } else {
        stop("Dimension mismatch: x$aggr$quantiles must have ", n, " rows")
      }
    }
  }
  
  # Plot mean predictions
  plot(1:n, x$aggr$mean, type = "l", ylab = "Predicted Values", xlab = "Observation Index", 
       col = "black", lwd = 2, ...)
  
  # Plot quantiles with varying shades of grey
  if (!is.null(x$aggr$quantiles)) {
    n_quantiles <- ncol(x$aggr$quantiles)
    grey_shades <- grey.colors(n_quantiles, start = 0.8, end = 0.2)  # Light to dark grey
    matlines(1:n, x$aggr$quantiles, lty = 2, col = grey_shades)
    # Add legend for quantiles
    if (!is.null(colnames(x$aggr$quantiles))) {
      legend("topright", legend = colnames(x$aggr$quantiles), col = grey_shades, 
             lty = 2, cex = 0.8, bg = "white")
    } else {
      legend("topright", legend = paste0("Quantile ", seq_len(n_quantiles)), 
             col = grey_shades, lty = 2, cex = 0.8, bg = "white")
    }
  }
  invisible(NULL)
}
