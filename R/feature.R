# Title     : Features for use in GMJMCMC
# Objective : Define the features and get a useful way to turn them into a string
# Created by: jonlachmann
# Created on: 2021-02-10

# transform = 0 is a multiplication of the (two) features listed in f
# transform > 0 is a transform of the type denoted by the transform variable
# f is a list of features involved in the feature
# alpha is the coefficient before each feature listed in f,
# and also possibly one more for an intercept

# A feature matrix has the structure
# | TRANSFORM   WIDTH ALPHA0 |
# | DEPTH       FEAT1 ALPHA1 |
# | NA          FEAT2 ALPHA2 |
# | NA          FEAT3 ALPHA3 |
# ...

#' Create Method for \"feature\" Class
#'
#' @param transform A numeric denoting the transform type
#' @param features A list of features to include
#' @param trans.priors A vector of prior inclusion penalties for the different transformations.
#' @param alphas A numeric vector denoting the alphas to use
#' @noRd
create.feature <- function (transform, features, trans.priors, alphas=NULL) {
  # Given no alphas, assume no intercept and unit coefficients
  if (is.null(alphas)) alphas <- c(0, rep(1, length(features)))
  if (length(alphas) != (length(features) + 1)) stop("Invalid alpha/feature count")
  # Calculate the depth, operation count and width of the new feature
  if (transform == 0) {
    depth <- 1 + depth.feature(features[[1]]) + depth.feature(features[[2]])
    oc <- 1 + oc.feature(features[[1]]) + oc.feature(features[[2]])
    width <- 2 + width.feature(features[[1]]) + width.feature(features[[2]])
  }
  else {
    depth <- 0 # Assume 0 depth to find the deepest included feature
    oc <- length(features) - 1 # Every + is an operation, and the outer transform is also one
    oc <- oc + trans.priors[transform]
    width <- length(features) # Width is the number of features
    for (i in 1:(length(features))) {
      locdepth <- depth.feature(features[[i]])
      lococ <- oc.feature(features[[i]])
      locwidth <- width.feature(features[[i]])
      width <- width + locwidth
      oc <- oc + lococ
      if (locdepth > depth) depth <- locdepth
    }
    depth <- depth + 1 # Add 1 to depth to account for the outer transformation
  }

  # Generate the new feature matrix
  newFeature <- list(matrix(c(transform, depth, rep(NA,length(alphas)-2),
                      width, 1:(length(features)), alphas), length(alphas)))
  attr(newFeature[[1]], "oc") <- oc
  feature <- append(newFeature, features, 0)
  class(feature) <- "feature"
  return(feature)
}

#' Update Alphas on a Feature
#'
#' @param feature The feature to be updated
#' @param alphas The alphas that will be used
#' @param recurse If we are recursing, to note the number of alphas used
#' @noRd
update.alphas <- function (feature, alphas, recurse=FALSE) {
  feat <- feature[[length(feature)]]
  alpha <- 0
  # This is a more complex feature
  if (is.matrix(feat)) {
    # Adjust intercept if it is not multiplication
    if (feat[1,1] > 0 && nrow(feat) > 2) {
      alpha <- alpha + 1
      feat[1,3] <- alphas[alpha]
    }
    for (i in 2:nrow(feat)) {
      # Multiplication does not have alphas, and a zero intercept is no intercept
      if (feat[1,1] > 0 && feat[1,3] != 0) {
        alpha <- alpha + 1
        feat[i,3] <- alphas[alpha]
      }
      # If we have a nested feature, recurse into it
      if (is.list(feature[[feat[i,2]]])) {
        recur <- update.alphas(feature[[feat[i,2]]], alphas[alpha + seq_along(alphas)], TRUE)
        feature[[feat[i,2]]] <- recur$feature
        alpha <- alpha + recur$alpha
      }
    }
  }
  feature[[length(feature)]] <- feat
  if (recurse) return(list(alpha=alpha, feature=feature))
  else return(feature)
}

#' Print Method for \"feature\" Class
#'
#' @param x An object of class "feature"
#' @param dataset Set the regular covariates as columns in a dataset
#' @param fixed How many of the first columns in dataset are fixed and do not contribute to variable selection
#' @param alphas Print a "?" instead of actual alphas to prepare the output for alpha estimation
#' @param labels Should the covariates be named, or just referred to as their place in the data.frame.
#' @param round Should numbers be rounded when printing? Default is FALSE, otherwise it can be set to the number of decimal places.
#' @param ... Not used.
#' 
#' @return String representation of a feature
#' 
#' @examples
#' result <- gmjmcmc(x = matrix(rnorm(600), 100),
#' y = matrix(rnorm(100), 100), 
#' P = 2, 
#' transforms = c("p0", "exp_dbl"))
#' print(result$populations[[1]][1])
#' 
#' @export
print.feature <- function (x, dataset = FALSE, fixed = 0, alphas = FALSE, labels = FALSE, round = FALSE, ...) {
  fString <- ""
  feat <- x[[length(x)]]
  # This is a more complex feature
  if (is.matrix(feat)) {
    transforms <- getOption("gmjmcmc-transformations")
    if (is.null(transforms)) stop("Please set the gmjmcmc-transformations option to your non-linear functions (see ?set.transforms).")
    # Assume that we are not doing multiplication
    op <- "+"
    # Add the outer transform is there is one
    if (feat[1, 1] > 0) fString <- paste0(fString, transforms[[feat[1, 1]]], "(")
    # If g = 0, we are doing multiplication
    else {
      op <- "*"
      fString <- paste0(fString, "(")
    }
    # If we are printing rounded features for neat output, round all alphas
    if (round) {
      feat[, 3] <- round(feat[, 3], round)
    }
    for (j in seq_len(nrow(feat))) {
      # No plus or multiplication sign on the last one
      if (j == nrow(feat)) op <- ""
      # If this is an intercept just add it in
      if (j == 1 && feat[j, 3] != 0) {
        if (!alphas) fString <- paste0(fString, feat[j, 3], op)
        else fString <- paste0(fString, "?", op)
      }
      # Otherwise this is a feature or covariate, do a recursive conversion
      if (j != 1) {
        # Process alphas, which are only present if there is more than one term in the feature
        # this implies that the feature is not a multiplication (i.e. only one _term_).
        if ((nrow(feat) > 2 || feat[1, 3] != 0) && feat[1, 1] > 0) {
          if (alphas) fString <- paste0(fString, "?*")
          else fString <- paste0(fString, feat[j,3], "*")
        }
        fString <- paste0(fString, print.feature(x[[feat[j, 2]]], dataset, fixed, alphas, labels, round), op)
      }
    }
    fString <- paste0(fString, ")")
  }
  # This is a plain covariate
  else if (is.numeric(feat)) {
    if (dataset) fString <- paste0("data$x[,", fixed + feat, "]")
    else if (labels[1] != FALSE) fString <- labels[feat]
    else fString <- paste0("x", feat)
  } else stop("Invalid feature structure")
  return(fString)
}

# A function to get the depth of a feature
depth.feature <- function (feature) {
  feat <- feature[[length(feature)]]
  if (is.matrix(feat)) return(feat[2,1])
  else if (is.numeric(feat)) return(0)
  else stop("Invalid feature structure")
}

# A function to get the width of a feature
width.feature <- function (feature) {
  feat <- feature[[length(feature)]]
  if (is.matrix(feat)) return(feat[1,2])
  else if (is.numeric(feat)) return(1)
  else stop("Invalid feature structure")
}

# A function to get the oc (operation count) of a feature
oc.feature <- function (feature) {
  feat <- feature[[length(feature)]]
  if (is.matrix(feat)) return(attr(feat, "oc"))
  else if (is.numeric(feat)) return(0)
  else stop("Invalid feature structure")
}

# A function to get the complexity measures of a list of features
complex.features <- function (features) {
  featcount <- length(features)
  width <- rep(NA, featcount)
  oc <- rep(NA, featcount)
  depth <- rep(NA, featcount)
  for (i in 1:featcount) {
    width[i] <- width.feature(features[[i]])
    oc[i] <- oc.feature(features[[i]])
    depth[i] <- depth.feature(features[[i]])
  }
  return(list(width=width, oc=oc, depth=depth))
}
