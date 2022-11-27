#' Predict using a BGNLM model.
#' @param model The model to use.
#' @param x The new data to use for the prediction, a matrix where each row is an observation.
#' @param link The link function to use
#' @param populations The populations to use
#' @param quantiles The quantiles to calculate credible intervals for the posterior moddes (in model space).
predict.bgnlm <- function (model, x, link=function(x) x, populations="last", quantiles=c(0.025, 0.5, 0.975)) {
  if (populations=="last") pops.use <- length(model$populations)
  else if (populations=="all") pops.use <- seq_len(length(model$populations))
  else if (populations=="best") pops.use <- which.max(unlist(model$best.margs))

  # TODO: Generalize for multiple populations
  models <- model$models[[pops.use]]
  features <- model$populations[[pops.use]]

  # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
  x.precalc <- precalc.features(cbind(0, 1, x), features)[, -1]

  # Get the renormalized marginal probabilities and indices of non-duplicate models
  model.probs <- marginal.probs.renorm(models, "models")
  models <- models[model.probs$idx]

  yhat <- matrix(NA, nrow=nrow(x), ncol=length(models))
  for (i in seq_along(models)) {
    yhat[, i] <- link(x.precalc[, c(TRUE, models[[i]]$model), drop=FALSE] %*% models[[i]]$coefs)
  }

  mean.pred <- rowSums(yhat * as.numeric(model.probs$probs))
  quantiles <- apply(yhat, 1, weighted.quantiles, weights=model.probs$probs, prob=quantiles)

  return(list(mean=mean.pred, quantiles=quantiles))
}

#' Calculate weighted quantiles
#' @param values The values to use
#' @param weights The weights of the values
#' @param prob The probabilities of the quantiles to use
#' @return
weighted.quantiles <- function (values, weights, prob=c(0.025, 0.975)) {
  ordered <- order(values)
  P <- cumsum(weights[ordered])

  iv <- integer(length(prob))
  for (i in seq_along(iv)) {
    iv[i] <- which.max(P >= prob[i])
  }
  {values[ordered]}[iv]
}
