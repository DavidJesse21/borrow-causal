options(box.path = "R")

box::use(
  borrow/ps_mem[borrow_ps_mem],
  borrow/bart[borrow_bart],
  borrow/commensurate[borrow_commensurate],
  simfuns/utils[trycatch_log]
)

box::use(
  stats[lm, coef, confint, pt, pnorm, qnorm, gaussian, quantile, median, vcov],
  posterior[summarise_draws, rhat],
  data.table[...],
  withr[with_preserve_seed]
)


#' Borrowing algorithms
#' 
#' @param .data (`data.table()`)\cr
#'   The randomly generated data.
#' @param constants (`list()`)\cr
#'   A list of objects that are specified for all simulations.
#'   Currently, these include:
#'   * `ci_level`: The nominal level of the two-sided confidence/credible interval
#'   * `num_warmup`, `num_posterior`, `num_chains`: Specifications for MCMC sampling.
#' @param stan_model (`CmdStanModel`)\cr
#'   A Stan model object that has already been compiled.
#' @param stan_output_dir (`character(1)`)\cr
#'   Directory to which to save Stan output (posterior draws).
#'   See documentation of `$sample()` method for more details.
#' @param algo_ids (`numeric()`)\cr
#'   Used by `wrapper_all_algos()` only.
#'   An integer vector specifying the IDs of the algorithms that should be applied to the data.
#' 
#' @returns (`numeric(6)`)\cr
#'   A named numeric vector containing the (1) (posterior) point estimate of the treatment effect,
#'   (2) the corresponding (posterior) p-value, (3 and 4) the lower and upper bound of the 
#'   corresponding confidence/credible interval, (5) the Rhat convergence diagnostic, 
#'   and (6) a borrowing summary statistic whose meaning is specific to the respective method.
#'   Note that (5) and (6) might also contain NAs for some methods (not applicable).
#' 
#' @name algorithms


#' @rdname algorithms
wrapper_lm = function(.data, constants) {
  m = lm(y ~ trt, data = .data)
  
  est = coef(m)[["trt"]]
  se_est = sqrt(diag(vcov(m)))[["trt"]]
  test_stat = est / se_est
  # One-sided p-value for testing H0: b <= 0
  pval = pt(test_stat, df = m$df.residual, lower.tail = FALSE)
  # Two-sided confidence interval
  ci = unname(confint(m, level = constants$ci_level)["trt", ])
  
  out = c(est, pval, ci, NA_real_, NA_real_)
  names(out) = c("est", "pval", "ci_lower", "ci_upper", "rhat", "borrow_stat")
  
  return(out)
}


#' @rdname algorithms
#' @export
wrapper_current = function(.data, constants) {
  wrapper_lm(.data[source == 0], constants)
}


#' @rdname algorithms
#' @export
wrapper_pooled = wrapper_lm


#' @rdname algorithms
#' @export
wrapper_ps_mem = function(.data, constants) {
  m = borrow_ps_mem(
    data = .data,
    outcome = "y", treatment = "trt", source = "source",
    covariates = c("x1", "x2", "x3"),
    family = gaussian()
  )
  
  post_mean = m$posterior$effect[["mu"]]
  post_sd = sqrt(m$posterior$effect[["sigma2"]])
  
  pval = 1 - pnorm(0, mean = post_mean, sd = post_sd, lower.tail = FALSE)
  ci_alpha = (1 - constants$ci_level) / 2
  ci_probs = c(ci_alpha, 1 - ci_alpha)
  ci = qnorm(ci_probs, mean = post_mean, sd = post_sd)
  
  out = c(post_mean, pval, ci, NA_real_, m$mem$exchange_probs[["source1"]])
  names(out) = c("est", "pval", "ci_lower", "ci_upper", "rhat", "borrow_stat")
  
  return(out)
}


#' @rdname algorithms
#' @export
wrapper_bart = function(.data, constants) {
  m = with_preserve_seed({
    borrow_bart(
      data = .data,
      outcome = "y", treatment = "trt", source = "source",
      covariates = c("x1", "x2", "x3"),
      family = gaussian(),
      num_warmup = constants$num_warmup,
      num_posterior = constants$num_posterior,
      num_chains = constants$num_chains,
      seed = 123,
      verbose = FALSE,
      n.threads = 1
    )
  })
  
  ci_alpha = (1 - constants$ci_level) / 2
  ci_probs = c(ci_alpha, 1 - ci_alpha)
  
  out = summarise_draws(
    m$draws$ate,
    mean,
    \(x) 1 - mean(x > 0),
    \(x) quantile(x, probs = ci_probs),
    rhat
  )[1, 2:6] |>
    unlist(recursive = FALSE, use.names = FALSE)
  out = c(out, NA_real_)
  names(out) = c("est", "pval", "ci_lower", "ci_upper", "rhat", "borrow_stat")
  
  return(out)
}



#' @rdname algorithms
#' @export
wrapper_ps_commensurate = function(.data,
                                   constants,
                                   stan_model,
                                   stan_output_dir = NULL) {
  m = with_preserve_seed({
    borrow_commensurate(
      model = stan_model,
      data = .data,
      outcome = "y", treatment = "trt", source = "source",
      covariates = c("x1", "x2", "x3"),
      num_warmup = constants$num_warmup,
      num_posterior = constants$num_posterior,
      num_chains = constants$num_chains,
      seed = 123,
      verbose = FALSE,
      parallel_chains = 1,
      output_dir = stan_output_dir
    )
  })
  
  ci_alpha = (1 - constants$ci_level) / 2
  ci_probs = c(ci_alpha, 1 - ci_alpha)
  
  out = summarise_draws(
    m$draws(variables = "beta_trt"),
    mean,
    \(x) 1 - mean(x > 0),
    \(x) quantile(x, probs = ci_probs),
    rhat
  )[1, 2:6] |>
    unlist(recursive = FALSE, use.names = FALSE)
  
  comm_param_post_median = summarise_draws(
    m$draws(variables = "tau"),
    median
  )[1, 2] |>
    unlist(recursive = FALSE, use.names = FALSE)
  
  out = c(out, comm_param_post_median)
  names(out) = c("est", "pval", "ci_lower", "ci_upper", "rhat", "borrow_stat")
  
  return(out)
}


#' @rdname algorithms
#' 
#' @note
#' A list helper object for iterating through the different algorithms
#' 
#' @export
li_algos = list(
  current = wrapper_current,
  pooled = wrapper_pooled,
  ps_mem = wrapper_ps_mem,
  ps_commensurate = wrapper_ps_commensurate,
  bart = wrapper_bart
)

#' @rdname algorithms
#' 
#' @note
#' A data.table helper object mapping IDs to the names of the algorithms.
#' 
#' @export
dt_algos = data.table(
  algo.id = 1:5,
  algo.name = c("current", "pooled", "ps_mem", "ps_commensurate", "bart")
)


#' @rdname algorithms
#' @export
wrapper_all_algos = function(.data,
                             constants,
                             stan_model,
                             stan_output_dir,
                             algo_ids = 1:5) {
  # Get all selected algorithms/functions
  algo_ids = sort(algo_ids)
  algo_names = dt_algos[algo.id %in% algo_ids, algo.name]
  algos = li_algos[algo_names]
  
  # Safely apply each algorithm to the generated data set
  res = lapply(names(algos), function(algo) {
    .fun = algos[[algo]]
    if (algo == "ps_commensurate") {
      trycatch_log(.fun(.data, constants, stan_model, stan_output_dir), err_val = rep(NA_real_, 6))
    } else {
      trycatch_log(.fun(.data, constants), err_val = rep(NA_real_, 6))
    }
  })
  
  # Combine, organize, and return results
  res = as.data.table(do.call(rbind, res))
  res[, algo.id := as.integer(algo_ids)]
  return(res)
}
