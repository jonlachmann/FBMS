# Title     : MJMCMC
# Objective : Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-04-27

#' Main algorithm for MJMCMC
#'
#' @param data A matrix containing the data to use in the algorithm,
#' first column should be the dependent variable, second should be the intercept
#' and the rest of the columns should be the independent variables.
#' @param loglik.pi The (log) density to explore
#' @param N The number of iterations to run for
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#'
#' @export mjmcmc
mjmcmc <- function (data, loglik.pi, N, probs, params, sub=F) {
  # Verify that data is well-formed
  data <- check.data(data)
  # Acceptance probability
  accept <- 0

  # Create a population of just the covariates
  S <- gen.covariates(ncol(data)-2)
  complex <- complex.features(S)

  # Initialize first model
  model.cur <- as.logical(rbinom(n = length(S), size = 1, prob = 0.5))
  model.cur.res <- loglik.pre(loglik.pi, model.cur, complex, data, params$loglik)
  model.cur <- list(prob=0, model=model.cur, coefs=model.cur.res$coefs, crit=model.cur.res$crit, alpha=0)
  best.crit <- model.cur$crit # Set first best criteria value

  cat("\nMJMCMC begin.\n")
  result <- mjmcmc.loop(data, complex, loglik.pi, model.cur, N, probs, params, sub)
  cat("\nMJMCMC done.\n")
  # Calculate acceptance rate
  result$accept <- result$accept / N
  result$populations <- S
  # Return formatted results
  return(result)
}

#' The main loop for the MJMCMC algorithm, used in both MJMCMC and GMJMCMC
#'
#' @param data The data to use
#' @param complex The complexity measures of the data
#' @param loglik.pi The (log) density to explore
#' @param model.cur The model to start the loop at
#' @param N The number of iterations to run for
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#'
#' @return A list containing the visited models, the models visited during local optimisation,
#' the acceptance count and the best critical value encountered.
mjmcmc.loop <- function (data, complex, loglik.pi, model.cur, N, probs, params, sub=F) {
  # Acceptance count
  accept <- 0
  # Number of covariates or features
  covar_count <- ncol(data) - 2
  # A list of models that have been visited
  models <- vector("list", N)
  # Initialize a vector to contain local opt visited models
  lo.models <- vector("list", 0)
  # Initialize list for keeping track of unique visited models
  visited.models <- list(models=matrix(model.cur$model, 1, covar_count), crit=model.cur$crit, count=1)
  best.crit <- model.cur$crit # Set first best criteria value

  progress <- 0
  mcmc_total <- as.numeric(model.cur$model)
  for (i in 1:N) {
    if (N > 40 && i %% floor(N/40) == 0) progress <- print.progressbar(progress, 40)

    if (i > params$burn_in) pip_estimate <- mcmc_total/i
    else pip_estimate <- rep(1/covar_count, covar_count)

    proposal <- mjmcmc.prop(data, loglik.pi, model.cur, complex, pip_estimate, probs, params, visited.models)
    if (proposal$crit > best.crit) {
      best.crit <- proposal$crit
      cat(paste("\rNew best crit:", best.crit, "\n"))
    }

    # If we did a large jump and visited models to save
    if (!is.null(proposal$models)) {
      lo.models <- c(lo.models, proposal$models)
      # If we are doing subsampling and want to update best mliks
      if (sub) {
        for (mod in seq_along(proposal$models)) {
          # Check if we have seen this model before
          mod.idx <- vec_in_mat(visited.models$models[1:visited.models$count,,drop=F], proposal$models[[mod]]$model)
          if (mod.idx == 0) {
            # If we have not seen the model before, add it
            visited.models$count <- visited.models$count + 1
            visited.models$crit <- c(visited.models$crit, proposal$models[[mod]]$crit)
            visited.models$models <- rbind(visited.models$models, proposal$models[[mod]]$model)
          } # This is a model seen before, set the best of the values available
          else visited.models$crit[mod.idx] <- max(proposal$models[[mod]]$crit, visited.models$crit[mod.idx])
        }
      }
      proposal$models <- NULL
    }
    if (sub) {
      # Check if we have seen this model before
      mod.idx <- vec_in_mat(visited.models$models[1:visited.models$count,,drop=F], proposal$model)
      if (mod.idx == 0) {
        # If we have not seen the model before, add it
        visited.models$count <- visited.models$count + 1
        visited.models$crit <- c(visited.models$crit, proposal$crit)
        visited.models$models <- rbind(visited.models$models, proposal$model)
      } # This is a model seen before, set the best of the values available
      else visited.models$crit[mod.idx] <- max (proposal$crit, visited.models$crit[mod.idx])
    }

    if (log(runif(1)) <= proposal$alpha) {
      model.cur <- proposal
      accept <- accept + 1
    }
    mcmc_total <- mcmc_total + model.cur$model
    # Add the current model to the list of visited models
    models[[i]] <- model.cur
  }
  return(list(models=models, accept=accept, lo.models=lo.models, best.crit=best.crit))
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
#'
mjmcmc.prop <- function (data, loglik.pi, model.cur, complex, pip_estimate, probs, params, visited.models=NULL) {
  l <- runif(1)
  if (l < probs$large) {
    ### Large jump

    ### Select kernels to use for the large jump
    q.l <- sample.int(n = 4, size = 1, prob = probs$large.kern) # Select large jump kernel
    q.o <- sample.int(n = 2, size = 1, prob = probs$localopt) # Select optimizer function
    q.r <- sample.int(n = 4, size = 1, prob = probs$random) # Select randomization kernel

    # Generate and do large jump
    large.jump <- gen.proposal(model.cur$model, params$large, q.l, NULL, pip_estimate) # Get the large jump
    chi.0.star <- xor(model.cur$model, large.jump$swap) # Swap large jump indices

    # Optimize to find a mode
    localopt <- local.optim(chi.0.star, data, loglik.pi, !large.jump$swap, complex, q.o, params) # Do local optimization
    chi.k.star <- localopt$model


    # Randomize around the mode
    proposal <- gen.proposal(chi.k.star, params$random, q.r, !large.jump$swap, pip_estimate, prob=T)
    proposal$model <- xor(chi.k.star, proposal$swap)

    # Do a backwards large jump and add in the kernel used in local optim to use the same for backwars local optim.
    chi.0 <- xor(proposal$model, large.jump$swap)

    # Do a backwards local optimization
    localopt2 <- local.optim(chi.0, data, loglik.pi, !large.jump$swap, complex, q.o, params, kern=localopt$kern)
    chi.k <- localopt2$model
    # TODO: We could compare if chi.k is reached by optimising from gamma (model.cur) as "intended"

    ### Calculate acceptance probability
    # Set up the parameters that were used to generate the proposal
    prop.params <- list(neigh.min=params$random$min, neigh.max=params$random$max, neigh.size=proposal$S)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, chi.k, q.r, prop.params, pip_estimate) # Get probability of gamma given chi.k

    # Store models visited during local optimization
    proposal$models <- c(localopt$models, localopt2$models)
  } else {
    ### Regular MH step
    # Select MH kernel
    q.g <- sample.int(n = 6, size = 1, prob = probs$mh)
    # Generate the proposal
    proposal <- gen.proposal(model.cur$model, params$mh, q.g, NULL, pip_estimate, prob=T)
    proposal$model <- xor(proposal$swap, model.cur$model)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, model.cur$model, q.g, params$mh, pip_estimate)
  }
  # Calculate log likelihoods for the proposed model
  proposal.res <- loglik.pre(loglik.pi, proposal$model, complex, data, params$loglik)
  proposal$crit <- proposal.res$crit

  # TODO: Compare to a list of best mliks for all visited models,
  # TODO: update that list if our estimate is better, otherwise update our estimate.
  # TODO: Save all models visited by local optim, and update the best mliks if we see one during local optim.
  # If we are running with subsampling, check the list for a better mlik
  if (!is.null(visited.models)) {
    mod.idx <- vec_in_mat(visited.models$models[1:visited.models$count,,drop=F], proposal$model)
    if (mod.idx != 0) proposal$crit <- max(proposal$crit, visited.models$crit[mod.idx])
  }

  # Calculate acceptance probability for proposed model
  proposal$alpha <- min(0, (proposal$crit + model.cur$prob) - (model.cur$crit + proposal$prob))

  ### Format results and return them
  proposal$swap <- NULL; proposal$S <- NULL
  proposal$coefs <- proposal.res$coefs
  return(proposal)
}