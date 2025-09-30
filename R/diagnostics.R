# Title     : Diagnostic functions
# Objective : Functions to do diagnostics on a GMJMCMC run
# Created by: jonlachmann
# Created on: 2021-02-24


#' Plot Convergence Diagnostics for GMJMCMC or GMJMCMC Merged Results
#'
#' Plots the convergence of summary statistics (e.g., median, mean) of log posteriors or marginal likelihoods
#' over populations for a GMJMCMC or GMJMCMC merged result object, with confidence intervals.
#'
#' @param res Object of class \code{gmjmcmc} or \code{gmjmcmc_merged} containing results from a GMJMCMC run or merged runs.
#' @param FUN Function to compute summary statistics (e.g., \code{median}, \code{mean}). Default is \code{median}.
#' @param conf Numeric; confidence level for intervals (e.g., 0.95 for 95%). Default is 0.95.
#' @param burnin Integer; number of initial populations to skip. Default is 0.
#' @param window Integer; size of the sliding window for computing standard deviation. Default is 5.
#' @param ylim Numeric vector; y-axis limits for the plot. If \code{NULL}, computed from confidence intervals.
#' @param ... Additional graphical parameters passed to \code{plot} and \code{lines} (e.g., \code{col}, \code{lwd}, \code{lty}, \code{main}, \code{xlab}, \code{ylab}).
#'
#' @return Returns \code{invisible(NULL)}. The function is called for its side effect of producing a plot.
#'
#' @examples
#' data(exoplanet)
#' result <- fbms(semimajoraxis ~ ., data = exoplanet, method = "gmjmcmc", transforms = c("sin"))
#' diagn_plot(result, FUN = median, conf = 0.95, main = "Convergence Plot")
#'
#' @export
diagn_plot <- function(res, FUN = median, conf = 0.95, burnin = 0, window = 5, ylim = NULL, ...) {
  # Input validation
  stopifnot(
    "res must be of class 'gmjmcmc' or 'gmjmcmc_merged'" = inherits(res, c("gmjmcmc", "gmjmcmc_merged")),
    "res must contain best.log.posteriors or best.margs" = length(res$thread.best) > 0 || length(res$best.margs) > 0,
    "FUN must be a function" = is.function(FUN),
    "conf must be between 0 and 1" = is.numeric(conf) && conf > 0 && conf < 1,
    "burnin must be a non-negative integer" = is.numeric(burnin) && burnin >= 0 && burnin %% 1 == 0,
    "window must be a positive integer" = is.numeric(window) && window >= 1 && window %% 1 == 0
  )
  
  args <- list(...)
  args[["..."]] <- NULL  # Remove any "..." element to avoid warning
  
  # Extract matrix.results
  if (length(res$thread.best) > 0) {
    matrix.results <- res$best.log.posteriors
  } else {
    matrix.results <- as.matrix(unlist(res$best.margs))
  }
  stopifnot("matrix.results must have at least burnin + 1 rows" = nrow(matrix.results) > burnin)
  
  # Compute summary statistics and confidence intervals
  sr <- sapply((burnin + 1):nrow(matrix.results), function(x) FUN(matrix.results[x, ]))
  sds <- c(0, sapply(2:length(sr), function(x) sd(sr[max(1, x - window):x])))
  ub <- sr + qnorm(1 - (1 - conf) / 2) * sds
  lb <- sr - qnorm(1 - (1 - conf) / 2) * sds
  
  # Set default plot parameters
  if (is.null(ylim) && is.null(args$ylim)) {
    ylim <- c(min(lb), max(ub))
  } else if (is.null(ylim)) {
    ylim <- args$ylim
  }
  main <- if (!is.null(args$main)) args$main else "Convergence"
  xlab <- if (!is.null(args$xlab)) args$xlab else "Population"
  ylab <- if (!is.null(args$ylab)) args$ylab else "Summary"
  
  args$main <- NULL
  args$xlab <- NULL
  args$ylab <- NULL
  args$ylim <- NULL
  
  # Create plot
  do.call(plot, c(
    list(
      y = sr,
      x = (burnin + 1):nrow(matrix.results),
      type = "l",
      col = 1,
      ylim = ylim,
      main = main,
      xlab = xlab,
      ylab = ylab
    ),
    args
  ))
  do.call(lines, c(
    list(
      y = ub,
      x = (burnin + 1):nrow(matrix.results),
      col = 1,
      lty = 2
    ),
    args
  ))
  do.call(lines, c(
    list(
      y = lb,
      x = (burnin + 1):nrow(matrix.results),
      col = 1,
      lty = 2
    ),
    args
  ))
  
  invisible(NULL)
}


