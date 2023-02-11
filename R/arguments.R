# Title     : Arguments generator
# Objective : Functions intended to help the user generate proper arguments for the algorithm
# Created by: jonlachmann
# Created on: 2021-02-19

#' Generate a probability list for (G)MJMCMC
#'
#' @param transforms A list of the transformations used (to get the count),
#' default is false and gives a MJMCMC prob list.
#'
#' @export gen.probs.list
gen.probs.list <- function (transforms=F) {
  ### Create a probability list for algorithm

  ## Mode jumping algorithm probabilities
  large <- 0.05                         # probability of a large jump
  large.kern <- c(0, 0, 0, 1)           # probability for type of large jump, only allow type 1-4
  localopt.kern <- c(0.5, 0.5)          # probability for each localopt algorithm
  random.kern <- c(0.3, 0.3, 0.2, 0.2)  # probability for random jump kernels
  mh <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1) # probability for regular mh kernels

  # Compile the list
  probs <- list(large=large, large.kern=large.kern, localopt.kern=localopt.kern,
                random.kern=random.kern, mh=mh)

  ## Feature generation probabilities
  if (transforms[1] != F) {
    transcount <- length(transforms)
    filter <- 0.6                           # filtration threshold
    gen <- rep(1/4, 4)                      # probability for different feature generation methods
    trans <- rep(1/transcount, transcount)  # probability for each different nonlinear transformation

    probs$filter <- filter
    probs$gen <- gen
    probs$trans <- trans
  }

  return(probs)
}

#' Generate a parameter list for (G)MJMCMC
#'
#' @param data The dataset that will be used in the algorithm
#' @param G True if we want parameters for GMJMCMC, false for MJMCMC
#'
#' @export gen.params.list
gen.params.list <- function (data, G=F) {
  ### Create a list of parameters for the algorithm

  ## Get the dimensions of the data to set parameters based on it
  data.dim <- data.dims(data)
  ncov <- data.dim[2] - 2
  nobs <- data.dim[1]

  ## Local optimization parameters
  sa_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15),
                  neigh.size=1, neigh.min=1, neigh.max=2)               # Simulated annealing proposal kernel parameters
  sa_params <- list(t.init=10, t.min=0.0001, dt=3, M=12, kern=sa_kern)  # Simulated annealing parameters
  greedy_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15),
                      neigh.size=1, neigh.min=1, neigh.max=2)           # Greedy algorithm proposal kernel parameters
  greedy_params <- list(steps=20, tries=3, kern=greedy_kern)            # Greedy algorithm parameters (60 models default)

  ## MJMCMC parameters
  burn_in <- 100                                                        # TODO
  large_params <- list(neigh.size=as.integer(ncov*0.35),
                       neigh.min=as.integer(ncov*0.25),
                       neigh.max=as.integer(ncov*0.45))                 # Large jump parameters
  random_params <- list(neigh.size=1, neigh.min=1, neigh.max=2)         # Small random jump parameters
  mh_params <- list(neigh.size=1, neigh.min=1, neigh.max=2)             # Regular MH parameters
  ## Compile the list and return
  params <- list(burn_in=burn_in, mh=mh_params, large=large_params, random=random_params,
                 sa=sa_params, greedy=greedy_params, loglik=list())

  # Add GMJMCMC specific parameters
  if (G) {
    ## GM parameters
    feat_params <- list(D = 5, L = 15,                        # Hard limits on feature complexity
                        alpha = 0,                            # alpha strategy (0=None, 1,2,3=strategies as per Hubin et al.) TODO: Fully Bayesian
                        pop.max = as.integer(ncov * 1.5),     # Max features population size
                        keep.org = F,                         # Always keep original covariates in every population
                        prel.filter = 0,                      # Filtration threshold for first population (i.e. filter covariates even if keep.org=T)
                        keep.min = 0.8,                       # Minimum proportion of features to always keep [0,1]
                        eps = 0.05,                           # Inclusion probability limit for feature generation
                        check.col = T,                        # Whether the colinearity should be checked
                        max.proj.size = 15)                   # Maximum projection size
    params$feat <- feat_params
    params$rescale.large <- F

  }

  return(params)
}