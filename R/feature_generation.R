# Title     : Features generation for use in GMJMCMC
# Objective : Generate features for use in GMJMCMC at the population transition step
# Created by: jonlachmann
# Created on: 2021-02-10

# Generate a multiplication feature
gen.multiplication <- function (features, marg.probs) {
  # generate two features to be multiplied
  feats <- sample.int(n = length(features), size = 2, prob = marg.probs, replace = T)
  create.feature(0, features[feats])
}

# Generate a modification feature
gen.modification <- function (features, marg.probs, trans.probs) {
  feat <- sample.int(n = length(features), size = 1, prob = marg.probs)
  trans <- sample.int(n = length(trans.probs), size = 1, prob = trans.probs)
  create.feature(trans, features[feat])
}

# Generate a projection feature
gen.projection <- function (features, marg.probs, trans.probs, max.width, max.size) {
  if (!is.null(max.size)) {
    max.width <- min(max.width, max.size + 1)
  }
  feat.count <- sample.int(n = (min(max.width, (length(features)))-1), size = 1) # TODO: Should be a specific distribution?
  feats <- sample.int(n = length(features), size = feat.count, prob = marg.probs)
  trans <- sample.int(n = length(trans.probs), size = 1, prob = trans.probs)
  # TODO: Generate alphas properly using various methods
  alphas <- rep(1, length(feats)+1)
  create.feature(trans, features[feats], alphas)
}

# Generate a new features from the initial covariates
gen.new <- function (features, F.0.size) {
  covariate <- sample.int(n = F.0.size, size = 1)
  return(features[[covariate]])
}

# Select a feature to generate and generate it
gen.feature <- function (features, marg.probs, data, loglik.alpha, transforms, probs, F.0.size, params) {
  tries <- 0
  feat.ok <- F
  while (!feat.ok && tries < 50) {
    feat.type <- sample.int(n = 4, size = 1, prob = probs$gen)
    if (feat.type == 1) feat <- gen.multiplication(features, marg.probs)
    if (feat.type == 2) feat <- gen.modification(features, marg.probs, probs$trans)
    if (feat.type == 3) feat <- gen.projection(features, marg.probs, probs$trans, params$L, params$max.proj.size)
    if (feat.type == 4) feat <- gen.new(features, F.0.size)
    # Check that the feature is not too wide or deep
    if (!(depth.feature(feat) > params$D || width.feature(feat) > params$L)) {
      # Generate alphas using the strategy chosen
      if (params$alpha > 0) {
        feat <- gen.alphas(params$alpha, feat, data, loglik.alpha)
      }
      if (!is.null(feat)) {
        # Check for linear dependence of new the feature
        if (length(features) == F.0.size) feats <- list()
        else feats <- features[(F.0.size+1):length(features)]
        if (params$check.col && !check.collinearity(feat, feats, F.0.size))
          feat.ok <- T
        else if (!params$check.col)
          feat.ok <- T
      }
    }
    tries <- tries + 1
    params$eps <- min(params$eps+0.01, 0.5)
    marg.probs <- pmin(pmax(marg.probs, params$eps), (1-params$eps))
  }
  if (!feat.ok) return(NULL)
  else return(feat)
}

check.collinearity <- function (proposal, features, F.0.size) {
  # Add the proposal to the feature list for evaluation
  features[[length(features)+1]] <- proposal
  # Generate mock data to test with (avoiding too costly computations)
  mock.data <- matrix(c(runif((F.0.size*2), -100, 100), rep(1,F.0.size*2),
                        runif((F.0.size*2)*(F.0.size), -100, 100)), F.0.size*2, F.0.size+2)
  # Use the mock data to precalc the features
  mock.data.precalc <- precalc.features(mock.data, features)
  # Fit a linear model with the mock data precalculated features
  linearmod <- lm(as.data.frame(mock.data.precalc[,-2]))
  # Check if all coefficients were possible to calculate
  if (sum(is.na(linearmod$coefficients)) == 0) return(F)
  else return(T)
}

# Generate features to represent the covariates, just takes the count needed
gen.covariates <- function (count) {
  features <- list()
  for (i in 1:count) {
    features <- c(features, i)
    class(features[[i]]) <- "feature"
  }
  return(features)
}