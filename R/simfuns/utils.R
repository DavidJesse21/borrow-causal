#' Modified `tryCatch()` function
#' 
#' @description
#' In a long-running simulation study with multiple repetitions and multiple methods, 
#' it becomes likely that some instances (repetition, method) might fail.
#' However, this should not cause all other computations (queued or finished) to be discarded.
#' This function tries to achieve two goals:
#' 1) Capture errors and return informative values in all cases.
#' 2) Still print errors (and warnings) to the console so they get logged and can be 
#'    investigated later on.
#' 
#' @param code (`expression()`)\cr
#'   Code to be evaluated.
#' @param err_val (`Any`)\cr
#'   Return value in case of an error.
#'   
#' @returns (`Any`)\cr
#'   Either what the specified `code` returns or the specified `err_val` object.
#'   
#' @export
trycatch_log = function(code, err_val = NA_real_) {
  warn = NULL
  err = NULL
  
  value = withCallingHandlers(
    tryCatch(code, error = function(e) {
      err <<- e
      err_val
    }),
    warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    }
  )
  
  if (!is.null(warn)) cat_warning(warn)
  if (!is.null(err)) cat_error(err)
  
  return(value)
}

cat_error = function(err) {
  .call = deparse(err$call)
  .message = err$message
  txt = sprintf("Error in %s: %s\n", .call, .message)
  cat(txt)
}

cat_warning = function(warn) {
  .call = deparse(warn$call)
  .message = warn$message
  txt = sprintf("Warning in %s: %s\n", .call, .message)
  cat(txt)
}



#' Estimate and format the runtime of a job
#' 
#' @param secs_per_job (`numeric(1)`)\cr
#'   The (assumed/estimated) runtime for a single job in seconds
#' @param num_jobs (`numeric(1)`)\cr
#'   The number of replications of the job.
#' @param mult (`numeric(1)`)\cr
#'   A number that the runtime will be multiplied with for more conservative estimates.
#' @param parallel (`list()` or `NULL`)\cr
#'   Either `NULL` to estimate the sequential runtime or a named list with the elements:\cr
#'   * `num_cores` (`numeric(1)`): The number of available cores
#'   * `prop` (`numeric(1)`): The proportion of parallelizable jobs
#' 
#' @details
#' For the estimation of the runtime using parallel processes Amdahl's law is used.
#' (https://en.wikipedia.org/wiki/Amdahl%27s_law)
#'   
#' @export
estimate_runtime = function(secs_per_job,
                            num_jobs = 1000L,
                            mult = 1,
                            parallel = NULL) {
  total = secs_per_job * num_jobs * mult
  
  if (is.null(parallel)) {
    hours = floor(total / 3600)
    minutes = floor((total %% 3600) / 60)
    seconds = ceiling(total %% 60)
    fmt = sprintf("%02d:%02d:%02d", hours, minutes, seconds)
    return(fmt)
  } else {
    speedup = 1 / ((1 - parallel$prop) + (parallel$prop / parallel$num_cores))
    total = total / speedup
    hours = floor(total / 3600)
    minutes = floor((total %% 3600) / 60)
    seconds = ceiling(total %% 60)
    fmt = sprintf("%02d:%02d:%02d", hours, minutes, seconds)
    return(fmt)
  }
}
