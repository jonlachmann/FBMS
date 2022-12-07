#' Run multiple mjmcmc runs in parallel, merging the results before returning.
#' @param runs The number of runs to run
#' @param cores The number of cores to run on
#' @param ... Parameters to pass to mjmcmc
#' @return Merged results from multiple mjmcmc runs
mjmcmc.parallel <- function (runs, cores=getOption("mc.cores", 2L), ...) {
  results <- mclapply(seq_len(runs), function (x) { gmjmcmc(...) }, mc.cores=cores)
}


#' Run multiple gmjmcmc runs in parallel returning a list of all results.
#' @param runs The number of runs to run
#' @param cores The number of cores to run on
#' @param ... Parameters to pass to gmjmcmc
#' @return Results from multiple gmjmcmc runs
gmjmcmc.parallel <- function (runs, cores=getOption("mc.cores", 2L), ...) {
  results <- mclapply(seq_len(runs), function (x) { gmjmcmc(...) }, mc.cores=cores)
}
