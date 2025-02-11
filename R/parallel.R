#' rmclapply: Cross-platform mclapply/forking hack for Windows
#'
#' This function applies a function in parallel to a list or vector (`X`) using multiple cores.
#' On Linux/macOS, it uses `mclapply`, while on Windows it uses a hackish version of parallelism.
#' The Windows version is based on `parLapply` to mimic forking following Nathan VanHoudnos.
#' @param runs The runs to run
#' @param args The arguments to pass to fun
#' @param fun The function to run
#' @param mc.cores Number of cores to use for parallel processing. Defaults to `detectCores()`.
#'
#' @return A list of results, with one element for each element of `X`.
rmclapply <- function(runs, args, fun, mc.cores = NULL) {
  if (is.null(args$verbose)) args$verbose <- TRUE

  if (is.null(mc.cores)) {
    mc.cores <- min(length(runs), detectCores())
  }

  if (Sys.info()[["sysname"]] == "Windows" & mc.cores > 1) {
    if (args$verbose) {
      message("Using parallelization hack for Windows with parLapply.")
    }


    ## N.B. setting outfile to blank redirects output to
    ##      the master console, as is the default with
    ##      mclapply() on Linux / Mac
    cl <- makeCluster(mc.cores, outfile = "")
    loaded.package.names <- c(
      sessionInfo()$basePkgs,
      names(sessionInfo()$otherPkgs)
    )
    
    tryCatch({
      ## Copy over all of the objects within scope to
      ## all clusters.
      clusterEvalQ(cl, library(FBMS))
      clusterExport(cl, "args")
      clusterExport(cl, ls(all.names = TRUE, envir = globalenv()), envir = globalenv())
      # Load required packages on each cluster node
      parLapply(cl, seq_along(cl), function(xx) {
        lapply(loaded.package.names, function(pkg) {
          require(pkg, character.only = TRUE)
        })
      })
      ## Run the lapply in parallel
      res <- parLapply(cl, runs, function(x) {
        set.seed(NULL)
        set.seed(as.integer(x) + sample.int(100000,1))
        do.call(fun, args)
      })
      gc()
      return(res)
    }, finally = {
      ## Stop the cluster
      stopCluster(cl)
    })

    ## Warn the user if they are using Windows
    if (Sys.info()[["sysname"]] == "Windows" & args$verbose == TRUE) {
      message(paste(
        "\n",
        "   *** Microsoft Windows detected ***\n",
        "   \n",
        "   For technical reasons, the MS Windows version of mclapply()\n",
        "   is implemented as a serial function instead of a parallel\n",
        "   function.",
        "   \n\n",
        "   As a quick hack, we replace this serial version of mclapply()\n",
        "   with a wrapper to parLapply() for this R session. Please see\n\n",
        "     http://www.stat.cmu.edu/~nmv/2014/07/14/implementing-mclapply-on-windows \n\n",
        "   for details.\n\n"
      ))
    }
  } else {
    res <- mclapply(runs, function(x) do.call(fun, args), mc.cores = mc.cores)
    gc()
    return(res)
  }
}



#' Run multiple mjmcmc runs in parallel, merging the results before returning.
#' @param runs The number of runs to run
#' @param cores The number of cores to run on
#' @param ... Further parameters passed to mjmcmc.
#' @return Merged results from multiple mjmcmc runs
#'
#' @examples
#' result <- mjmcmc.parallel(runs = 1, cores = 1, matrix(rnorm(600), 100), gaussian.loglik)
#' summary(result)
#' plot(result)
#'
#' @export
mjmcmc.parallel <- function(runs = 2, cores = getOption("mc.cores", 2L), ...) {
  results <- rmclapply(seq_len(runs), args = list(...), mc.cores = cores, fun = mjmcmc)
  class(results) <- "mjmcmc_parallel"
  return(results)
}


#' Run multiple gmjmcmc (Genetically Modified MJMCMC) runs in parallel returning a list of all results.
#' @param runs The number of runs to run
#' @param cores The number of cores to run on
#' @param merge.options A list of options to pass to the [merge_results()] function run after the
#' @inheritParams gmjmcmc
#' @param ... Further parameters passed to mjmcmc.
#' @return Results from multiple gmjmcmc runs
#'
#' @examples
#' result <- gmjmcmc.parallel(
#'   runs = 1,
#'   cores = 1,
#'   list(populations = "best", complex.measure = 2, tol = 0.0000001),
#'   matrix(rnorm(600), 100),
#'   P = 2,
#'   gaussian.loglik,
#'   loglik.alpha = gaussian.loglik.alpha,
#'   c("p0", "exp_dbl")
#' )
#'
#' summary(result)
#'
#' plot(result)
#'
#' @export
gmjmcmc.parallel <- function(runs = 2, cores = getOption("mc.cores", 2L), merge.options = list(populations = "best", complex.measure = 2, tol = 0.0000001), data, loglik.pi = gaussian.loglik, loglik.alpha = gaussian.loglik.alpha, transforms, ...) {
  options("gmjmcmc-transformations" = transforms)

  results <- rmclapply(seq_len(runs), args = list(data = data, loglik.pi = loglik.pi, loglik.alpha = loglik.alpha, transforms = transforms, ...), mc.cores = cores, fun = gmjmcmc)

  class(results) <- "gmjmcmc_parallel"
  merged <- merge_results(results, merge.options$populations, merge.options$complex.measure, merge.options$tol, data = data)
  gc()
  return(merged)
}
