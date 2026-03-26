# Title     : Diagnostic functions
# Objective : Functions to do diagnostics on a GMJMCMC run
# Created by: jonlachmann
# Created on: 2021-02-24

# TODO: Inter-thread variance comparison of the marginal log posterior of the best found model.

#' Plot Convergence Diagnostics for GMJMCMC or GMJMCMC Merged Results
#'
#' Plots the convergence of summary statistics (e.g., median, mean) of log posteriors or marginal likelihoods
#' over populations for a GMJMCMC or GMJMCMC merged result object, with support for various plot types.
#'
#' @param res Object of class \code{gmjmcmc} or \code{gmjmcmc_merged} containing results from a GMJMCMC run or merged runs.
#' @param FUN Function to compute summary statistics (e.g., \code{median}, \code{mean}). Default is \code{median}.
#' @param conf Numeric; confidence level for intervals (e.g., 0.95 for 95%). Default is 0.95.
#' @param burnin Integer; number of initial populations to skip. Default is 0.
#' @param window Integer; size of the sliding window for computing standard deviation. Default is 5.
#' @param ylim Numeric vector; y-axis limits for the plot. If \code{NULL}, computed from the data being plotted.
#' @param type String; type of plot to produce. Options are:
#' \itemize{
#'   \item \code{"convergence"}: (Default) Summary statistic with confidence intervals based on sliding window SD.
#'   \item \code{"min-mean-max"}: Min, mean, and max log-posteriors per population across threads.
#'   \item \code{"all-threads"}: Each thread's summary statistic (defined by \code{FUN}) as a separate time series.
#'   \item \code{"total-mass"}: Total log-posterior (sum of log-posteriors of all unique models) per population.
#' }
#' @param per_thread Logical; if \code{TRUE} and \code{type = "total-mass"}, plots mass for each thread individually.
#' @param pool Logical; if \code{TRUE}, statistics are calculated across all models pooled from all threads for each population.
#' @param mass Logical; if \code{TRUE}, total statistics are calculated for each thread under "min-mean-max" settings.
#' @param ... Additional graphical parameters passed to \code{plot} and \code{lines} (e.g., \code{col}, \code{lwd}, \code{lty}, \code{main}, \code{xlab}, \code{ylab}).
#'
#' @return Returns a list of summary statistics used for the plot.
#'
#' @examples
#' data(exoplanet)
#' result <- fbms(semimajoraxis ~ ., data = exoplanet, method = "gmjmcmc", transforms = c("sin"))
#' diagn_plot(result, FUN = median, conf = 0.95, main = "Convergence Plot")
#'
#' @export
diagn_plot <- function(res, FUN = median, conf = 0.95, burnin = 0, window = 5, ylim = NULL, type = "convergence", per_thread = FALSE, pool = FALSE, mass = FALSE, ...) {
  # Input validation
  stopifnot(
    "res must be of class 'gmjmcmc' or 'gmjmcmc_merged'" = inherits(res, c("gmjmcmc", "gmjmcmc_merged")),
    "FUN must be a function" = is.function(FUN),
    "conf must be between 0 and 1" = is.numeric(conf) && conf > 0 && conf < 1,
    "burnin must be a non-negative integer" = is.numeric(burnin) && burnin >= 0 && burnin %% 1 == 0,
    "window must be a positive integer" = is.numeric(window) && window >= 1 && window %% 1 == 0,
    "type must be one of 'convergence', 'min-mean-max', 'all-threads', 'total-mass'" = type %in% c("convergence", "min-mean-max", "all-threads", "total-mass")
  )
  
  args <- list(...)
  args[["..."]] <- NULL # Remove any "..." element to avoid warning
  
  # Extract results list
  if (inherits(res, "gmjmcmc_merged")) {
    results_list <- res$results.raw
  } else {
    results_list <- list(res)
  }
  
  num_pops_total <- length(results_list[[1]]$models)
  num_threads <- length(results_list)
  pops <- (burnin + 1):num_pops_total
  
  if (num_threads == 0) stop("No results available in the object.")
  
  # Aggregated statistics and values
  stat.matrix <- NULL
  pooled_mins <- NULL
  pooled_means <- NULL
  pooled_maxs <- NULL
  pooled_sr <- NULL
  pooled_masses <- NULL
  
  if (pool) {
    pooled_mins <- numeric(num_pops_total)
    pooled_means <- numeric(num_pops_total)
    pooled_maxs <- numeric(num_pops_total)
    pooled_sr <- numeric(num_pops_total)
    pooled_masses <- numeric(num_pops_total)
    
    for (p in 1:num_pops_total) {
      all_crits <- c()
      all_unique_models <- list()
      all_unique_crits <- c()
      
      for (t in 1:num_threads) {
        if (length(results_list[[t]]$models) >= p) {
          current_crits <- sapply(results_list[[t]]$models[[p]], function(x) x$crit)
          all_crits <- c(all_crits, current_crits)
          
          
          
          if (type == "total-mass") {
            idx <- results_list[[t]]$model.probs.idx[[p]]
            all_unique_models <- c(all_unique_models, results_list[[t]]$models[[p]][idx])
            all_unique_crits <- c(all_unique_crits, sapply(results_list[[t]]$models[[p]][idx], function(x) x$crit))
          }
        }
        if(p == 1 & t==1)
          maxcrit.1 <- max(all_unique_crits)
      }
      
      if (length(all_crits) > 0) {
        pooled_mins[p] <- min(all_crits)
        pooled_means[p] <- mean(all_crits)
        pooled_maxs[p] <- max(all_crits)
        pooled_sr[p] <- FUN(all_crits)
        
        if (type == "total-mass") {
          model_size <- length(all_unique_models[[1]]$model)
          model_matrix <- matrix(unlist(lapply(all_unique_models, function(x) x$model)), ncol = model_size, byrow = TRUE)
          duplicates <- duplicated(model_matrix)
          pooled_masses[p] <- sum(exp(all_unique_crits[!duplicates] - res$crit.best))
          
        }
      }
    }
  } else {
    # Per-thread statistics
    stat.matrix <- matrix(NA, nrow = num_pops_total, ncol = num_threads)
    for (t in 1:num_threads) {
      if (is.null(results_list[[t]]$models)) next
      for (p in 1:num_pops_total) {
        if (length(results_list[[t]]$models) >= p) {
          if (type == "total-mass") {
            unique_idx <- results_list[[t]]$model.probs.idx[[p]]
            crits <- sapply(results_list[[t]]$models[[p]][unique_idx], function(x) x$crit)
            stat.matrix[p, t] <- sum(exp(crits-res$crit.best))
          } else {
            unique_idx <- results_list[[t]]$model.probs.idx[[p]]
            crits <- sapply(results_list[[t]]$models[[p]][unique_idx], function(x) x$crit)
            if(mass)
            {
              if(p == 1 & t == 1)
              {  
                FUN = function(x){sum(exp(x-res$crit.best))}
              }
            }
              
            stat.matrix[p, t] <- FUN(crits)
          }
        }
      }
    }
  }
  
  # Default labels and titles
  main <- if (!is.null(args$main)) args$main else paste("Convergence (", type, if(pool) " pooled" else "", ")", sep = "")
  xlab <- if (!is.null(args$xlab)) args$xlab else "Population"
  ylab <- if (!is.null(args$ylab)) args$ylab else if (type == "total-mass" || mass) "Total mass (unique)" else "Log Posterior"
  
  # Helper to remove used plot args to avoid warnings in do.call
  clean_args <- function(a) {
    a$main <- NULL; a$xlab <- NULL; a$ylab <- NULL; a$ylim <- NULL; a$col <- NULL; a$lty <- NULL; a$type <- NULL
    return(a)
  }
  
  if (type == "convergence") {
    if (pool) {
      sr <- pooled_sr[pops]
    } else {
      sr <- apply(stat.matrix[pops, , drop = FALSE], 1, FUN)
    }
    
    sds <- c(0, sapply(2:length(sr), function(x) sd(sr[max(1, x - window):x])))
    if (length(sds) > 1 && sds[1] == 0) sds[1] <- sds[2]
    
    ub <- sr + qnorm(1 - (1 - conf) / 2) * sds
    lb <- sr - qnorm(1 - (1 - conf) / 2) * sds
    
    if (is.null(ylim)) ylim <- c(min(lb, na.rm = TRUE), max(ub, na.rm = TRUE))
    do.call(plot, c(list(y = sr, x = pops, type = "l", col = 1, ylim = ylim, main = main, xlab = xlab, ylab = ylab), clean_args(args)))
    do.call(lines, c(list(y = ub, x = pops, col = 1, lty = 2), clean_args(args)))
    do.call(lines, c(list(y = lb, x = pops, col = 1, lty = 2), clean_args(args)))
    return(invisible(list(stat = sr, lower = lb, upper = ub)))
    
  } else if (type == "min-mean-max") {
    if (pool) {
      mins <- pooled_mins[pops]
      means <- pooled_means[pops]
      maxs <- pooled_maxs[pops]
    } else {
      mins <- apply(stat.matrix[pops, , drop = FALSE], 1, min)
      means <- apply(stat.matrix[pops, , drop = FALSE], 1, mean)
      maxs <- apply(stat.matrix[pops, , drop = FALSE], 1, max)
    }
    
    if (is.null(ylim)) ylim <- c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))
    do.call(plot, c(list(y = means, x = pops, type = "l", col = 1, ylim = ylim, main = main, xlab = xlab, ylab = ylab), clean_args(args)))
    do.call(lines, c(list(y = mins, x = pops, col = "blue", lty = 2), clean_args(args)))
    do.call(lines, c(list(y = maxs, x = pops, col = "blue", lty = 2), clean_args(args)))
    legend("bottomright", legend = c("Max", "Mean", "Min"), col = c("blue", 1, "blue"), lty = c(2, 1, 2))
    return(invisible(list(min = mins, mean = means, max = maxs)))
    
  } else if (type == "all-threads") {
    if (pool) stop("all-threads cannot be used with pool = TRUE as it collapses threads.")
    mat <- stat.matrix[pops, , drop = FALSE]
    if (is.null(ylim)) ylim <- c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE))
    do.call(plot, c(list(y = mat[, 1], x = pops, type = "l", col = 1, ylim = ylim, main = main, xlab = xlab, ylab = ylab), clean_args(args)))
    if (ncol(mat) > 1) {
      for (i in 2:ncol(mat)) do.call(lines, c(list(y = mat[, i], x = pops, col = i), clean_args(args)))
    }
    return(invisible(mat))
    
  } else if (type == "total-mass") {
    if (pool) {
      masses <- pooled_masses[pops]
      if (is.null(ylim)) ylim <- c(min(masses, na.rm = TRUE), max(masses, na.rm = TRUE))
      do.call(plot, c(list(y = masses, x = pops, type = "l", col = 4, ylim = ylim, main = main, xlab = xlab, ylab = ylab), clean_args(args)))
      return(invisible(list(mass = masses)))
    } else {
      mat <- stat.matrix[pops, , drop = FALSE]
      if (per_thread) {
        if (is.null(ylim)) ylim <- c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE))
        do.call(plot, c(list(y = mat[, 1], x = pops, type = "l", col = 1, ylim = ylim, main = main, xlab = xlab, ylab = ylab), clean_args(args)))
        if (ncol(mat) > 1) {
          for (i in 2:ncol(mat)) do.call(lines, c(list(y = mat[, i], x = pops, col = i), clean_args(args)))
        }
        return(invisible(mat))
      } else {
        total_mass <- rowSums(mat, na.rm = TRUE)
        if (is.null(ylim)) ylim <- c(min(total_mass, na.rm = TRUE), max(total_mass, na.rm = TRUE))
        do.call(plot, c(list(y = total_mass, x = pops, type = "l", col = 4, ylim = ylim, main = main, xlab = xlab, ylab = ylab), clean_args(args)))
        return(invisible(list(mass = total_mass)))
      }
    }
  }
}
