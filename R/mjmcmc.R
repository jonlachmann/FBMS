# Title     : MJMCMC
# Objective : Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-04-27

#' Main algorithm for MJMCMC (Genetically Modified MJMCMC)
#'
#' @param data A matrix containing the data to use in the algorithm,
#' first column should be the dependent variable, second should be the intercept
#' and the rest of the columns should be the independent variables.
#' @param loglik.pi The (log) density to explore
#' @param N The number of iterations to run for
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#' @param verbose A logical denoting if messages should be printed
#'
#' @return A list containing the following elements:
#' \item{models}{All visited models.}
#' \item{accept}{Average acceptance rate of the chain.}
#' \item{lo.models}{All models visited during local optimization.}
#' \item{best.crit}{The highest log marginal probability of the visited models.}
#' \item{marg.probs}{Marginal probabilities of the features.}
#' \item{model.probs}{Marginal probabilities of all of the visited models.}
#' \item{model.probs.idx}{Indices of unique visited models.}
#' \item{populations}{The covariates represented as a list of features.}
#'
#' @examples
#' result <- mjmcmc(matrix(rnorm(600), 100), gaussian.loglik)
#' summary(result)
#' plot(result)
#'
#' @export mjmcmc
mjmcmc <- function (data, loglik.pi = gaussian.loglik, N = 100, probs = NULL, params = NULL, sub = FALSE, verbose = TRUE) {
  # Verify that data is well-formed
  data <- check.data(data, verbose)

  # Generate default probabilities and parameters if there are none supplied.
  if (is.null(probs)) probs <- gen.probs.mjmcmc()
  if (is.null(params)) params <- gen.params.mjmcmc(data)

  # Acceptance probability
  accept <- 0

  # Create a population of just the covariates
  S <- gen.covariates(ncol(data)-2)
  complex <- complex.features(S)

  # Initialize first model
  model.cur <- as.logical(rbinom(n = length(S), size = 1, prob = 0.5))
  model.cur.res <- loglik.pre(loglik.pi, model.cur, complex, data, params$loglik, visited.models=NULL, sub = sub)
  model.cur <- list(prob=0, model=model.cur, coefs=model.cur.res$coefs, crit=model.cur.res$crit, alpha=0)

  if (verbose) cat("\nMJMCMC begin.\n")
  result <- mjmcmc.loop(data, complex, loglik.pi, model.cur, N, probs, params, sub, verbose)
  if (verbose) cat("\nMJMCMC done.\n")
  # Calculate acceptance rate
  result$accept <- result$accept / N
  result$populations <- S

  # Return formatted results
  class(result) <- "mjmcmc"
  return(result)
}

#' The main loop for the MJMCMC (Mode Jumping MCMC) algorithm, used in both MJMCMC and GMJMCMC (Genetically Modified MJMCMC)
#'
#' @param data The data to use
#' @param complex The complexity measures of the data
#' @param loglik.pi The (log) density to explore
#' @param model.cur The model to start the loop at
#' @param N The number of iterations to run for
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#' @param verbose A logical denoting if messages should be printed
#'
#' @return A list containing the following elements:
#' \item{models}{All visited models.}
#' \item{accept}{Number of accepted proposals of the chain.}
#' \item{lo.models}{All models visited during local optimization.}
#' \item{best.crit}{The highest log marginal probability of the visited models.}
#' \item{marg.probs}{Marginal probabilities of the features.}
#' \item{model.probs}{Marginal probabilities of all of the visited models.}
#' \item{model.probs.idx}{Indices of unique visited models.}
#'
#' @noRd
#'
mjmcmc.loop <- function (data, complex, loglik.pi, model.cur, N, probs, params, sub = FALSE, verbose = TRUE) {
  # Acceptance count
  accept <- 0
  # Number of covariates or features
  covar_count <- ncol(data) - 2
  # A list of models that have been visited
  models <- vector("list", N)
  # Initialize a vector to contain local opt visited models
  lo.models <- vector("list", 0)
  # Initialize list for keeping track of unique visited models
  visited.models <- hashmap()
  visited.models[[model.cur$model]] <- list(crit = model.cur$crit, coefs = model.cur$coefs)
  best.crit <- model.cur$crit # Set first best criteria value

  progress <- 0
  mcmc_total <- as.numeric(model.cur$model)
  for (i in seq_len(N)) {
    if (verbose && N > 40 && i %% floor(N / 40) == 0) progress <- print_progressbar(progress, 40)

    if (i > params$burn_in) pip_estimate <- mcmc_total / i
    else pip_estimate <- rep(1 / covar_count, covar_count)

    proposal <- mjmcmc.prop(data, loglik.pi, model.cur, complex, pip_estimate, probs, params, visited.models, sub = sub)
    if (proposal$crit > best.crit) {
      best.crit <- proposal$crit
      if (verbose) cat(paste("\rNew best population crit:", best.crit, "\n"))
    }

    # If we did a large jump and visited models to save
    if (!is.null(proposal$models)) {
      lo.models <- c(lo.models, proposal$models)
      for (mod in seq_along(proposal$models)) {
        visited.models[[proposal$models[[mod]]$model]] <- list(crit = proposal$models[[mod]]$crit, coefs = proposal$models[[mod]]$coefs)
      }
      proposal$models <- NULL
    }
    visited.models[[proposal$model]] <- list(crit = proposal$crit, coefs = proposal$coefs)

    if (log(runif(1)) <= proposal$alpha) {
      model.cur <- proposal
      accept <- accept + 1
    }
    mcmc_total <- mcmc_total + model.cur$model
    # Add the current model to the list of visited models
    models[[i]] <- model.cur
  }

  # Calculate and store the marginal inclusion probabilities and the model probabilities
  marg.probs <- marginal.probs.renorm(c(models, lo.models), type = "both")

  return(list(
    models = models,
    accept = accept,
    lo.models = lo.models,
    best.crit = best.crit,
    marg.probs = marg.probs$probs.f,
    model.probs = marg.probs$probs.m,
    model.probs.idx = marg.probs$idx
  ))
}

#' Subalgorithm for generating a proposal and acceptance probability in (G)MJMCMC
#'
#' @param data The data to use in the algorithm
#' @param loglik.pi The the (log) density to explore
#' @param model.cur The current model to make the proposal respective to
#' @param complex The complexity measures used when evaluating the marginal likelihood
#' @param pip_estimate The current posterior inclusion probability estimate, used for proposals
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param visited.models A list of the previously visited models to use when subsampling and avoiding recalculation
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#'
#' @noRd
#'
mjmcmc.prop <- function (data, loglik.pi, model.cur, complex, pip_estimate, probs, params, visited.models = NULL, sub = FALSE) {
  l <- runif(1)
  if (l < probs$large) {
    ### Large jump

    ### Select kernels to use for the large jump
    q.l <- sample.int(n = 4, size = 1, prob = probs$large.kern) # Select large jump kernel
    q.o <- sample.int(n = 2, size = 1, prob = probs$localopt) # Select optimizer function
    q.r <- sample.int(n = 2, size = 1, prob = probs$random.kern) # Set randomization kernel

    # Generate and do large jump
    large.jump <- gen.proposal(model.cur$model, params$large, q.l, NULL, pip_estimate) # Get the large jump
    chi.0.star <- xor(model.cur$model, large.jump$swap) # Swap large jump indices

    # Optimize to find a mode
    localopt <- local.optim(chi.0.star, data, loglik.pi, !large.jump$swap, complex, q.o, params, visited.models = visited.models, sub = sub) # Do local optimization
    chi.k.star <- localopt$model

    # Randomize around the mode
    proposal <- gen.proposal(chi.k.star, list(neigh.size = length(pip_estimate), neigh.min = 1, neigh.max = length(pip_estimate)), q.r, NULL, (pip_estimate * 0 + 1 - params$random$prob), prob=TRUE)
    proposal$model <- xor(chi.k.star, proposal$swap)

    # Do a backwards large jump and add in the kernel used in local optim to use the same for backwards local optim.
    chi.0 <- xor(proposal$model, large.jump$swap)

    # Do a backwards local optimization
    localopt2 <- local.optim(chi.0, data, loglik.pi, !large.jump$swap, complex, q.o, params, kernel = localopt$kern,  visited.models=visited.models, sub = sub)
    chi.k <- localopt2$model

    ### Calculate acceptance probability
    # Set up the parameters that were used to generate the proposal
    prop.params <- list(neigh.size = length(pip_estimate), neigh.min = 1, neigh.max = length(pip_estimate))#list(neigh.min = params$random$min, neigh.max = params$random$max, neigh.size = proposal$S)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, chi.k, q.r, prop.params, (pip_estimate*0 + 1 - params$random$prob)) # Get probability of gamma given chi.k

    # Store models visited during local optimization
    proposal$models <- c(localopt$models, localopt2$models)
  } else {
    ### Regular MH step
    # Select MH kernel
    q.g <- sample.int(n = 6, size = 1, prob = probs$mh)
    # Generate the proposal
    proposal <- gen.proposal(model.cur$model, params$mh, q.g, NULL, pip_estimate, prob = TRUE)
    proposal$model <- xor(proposal$swap, model.cur$model)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, model.cur$model, q.g, params$mh, pip_estimate)
  }
  # Calculate log likelihoods for the proposed model
  proposal.res <- loglik.pre(loglik.pi, proposal$model, complex, data, params$loglik, visited.models=visited.models, sub = sub)
  proposal$crit <- proposal.res$crit

  # Calculate acceptance probability for proposed model
  proposal$alpha <- min(0, (proposal$crit + model.cur$prob) - (model.cur$crit + proposal$prob))

  ### Format results and return them
  proposal$swap <- NULL; proposal$S <- NULL
  proposal$coefs <- proposal.res$coefs
  return(proposal)
}
