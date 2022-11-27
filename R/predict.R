#' Predict using a BGNLM model.
#' @param model The model to use.
#' @param x The new data to use for the prediction, a matrix where each row is an observation.
#' @param link The link function to use
#' @param populations The populations to use
predict.bgnlm <- function (model, x, link=function(x) x, populations="last") {
  if (populations=="last") pops.use <- length(model$populations)
  else if (populations=="all") pops.use <- seq_len(length(model$populations))
  else if (populations=="best") pops.use <- which.max(unlist(model$best.margs))

   # TODO: Generalize for multiple populations
  models <- model$models[[pops.use]]
  features <- model$populations[[pops.use]]

  # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
  x.precalc <- precalc.features(cbind(0, 1, x), features)[, -1]

}

x <- kmdata[1:10, -1]
