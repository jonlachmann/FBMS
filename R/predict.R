#' Predict using a BGNLM model.
#'
#' @param model The model to use.
#' @param x The new data to use for the prediction, a matrix where each row is an observation.
#' @param link The link function to use
#' @param quantiles The quantiles to calculate credible intervals for the posterior moddes (in model space).
#'
#' @export predict.bgnlm
predict.bgnlm <- function (model, x, link=function(x) x, quantiles=c(0.025, 0.5, 0.975)) {
  x <- as.matrix(x)
  preds <- list()
  for (i in seq_along(model$results)) {
    preds[[i]] <- list()
    for (j in seq_along(model$results[[i]]$populations)) {
      # Select the models and features to predict from at this iteration
      models <- model$results[[i]]$models[[j]]
      features <- model$results[[i]]$populations[[j]]
      model.probs <- model$results[[i]]$model.probs[[j]]

      # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
      x.precalc <- precalc.features(cbind(0, 1, x), features)[, -1]

      yhat <- matrix(NA, nrow=nrow(x), ncol=length(models))
      for (k in seq_along(models)) {
        yhat[, k] <- link(x.precalc[, c(TRUE, models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
      }

      mean.pred <- rowSums(yhat %*% diag(as.numeric(model.probs)))
      pred.quant <- apply(yhat, 1, weighted.quantiles, weights=model.probs, prob=quantiles)

      preds[[i]][[j]] <- list(mean=mean.pred, quantiles=pred.quant, weight=model$results[[i]]$pop.weights[j])
    }
  }

  aggr <- list()
  aggr$mean <- 0 * preds[[1]][[1]]$mean
  aggr$quantiles <- 0 * preds[[1]][[1]]$quantiles
  for (i in seq_along(preds)) {
    for (j in seq_along(preds[[i]])) {
      aggr$mean <- aggr$mean + preds[[i]][[j]]$mean * model$results[[i]]$pop.weights[j]
      aggr$quantiles <- aggr$quantiles + preds[[i]][[j]]$quantiles * model$results[[i]]$pop.weights[j]
    }
  }

  return(list(aggr=aggr, preds=preds))
}

#' Calculate weighted quantiles
#'
#' @param values The values to use
#' @param weights The weights of the values
#' @param prob The probabilities of the quantiles to use
#'
#' @return Weighted quantiles
weighted.quantiles <- function (values, weights, prob=c(0.025, 0.975)) {
  ordered <- order(values)
  P <- cumsum(weights[ordered])

  iv <- integer(length(prob))
  for (i in seq_along(iv)) {
    iv[i] <- which.max(P >= prob[i])
  }
  {values[ordered]}[iv]
}
