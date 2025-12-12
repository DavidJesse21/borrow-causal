box::use(
  data.table[...],
  dbarts[bart2],
  posterior[as_draws],
  stats[binomial, gaussian, as.formula, pnorm, rgamma],
  chk = checkmate
)

box::use(
  ./utils[assert_binary]
)


#' Bayesian dynamic borrowing using BART
#' 
#' @description 
#' This function uses the BART model to borrow historical control data for treatment effect estimation based on the 
#' paper by Zhou and Ji (see references below).
#' Based on the potential outcomes framework and employing the Bayesian bootstrap it returns the 
#' (population) average treatment effect.
#' 
#' @param data (`data.frame()`)\cr
#'   A data.frame containing the combined current and historical trial data.
#' @param outcome (`character(1)`)\cr
#'   A string specifying the column name for the outcome variable.
#' @param treatment (`character(1)`)\cr
#'   A string specifying the column name for the treatment variable.
#' @param source (`character(1)`)\cr
#'   A string specifying the column name for the source/study indicator variable.
#' @param covariates (`character()`)\cr
#'   A character vector specifying the column names of covariates to be used in the models.
#' @param family (`family`)\cr
#'   A family object (see `?stats::family`).
#'   Possible options are `binomial()` with probit link and `gaussian()` with identity link.
#' @param num_warmup,num_posterior,num_chains,seed (`numeric(1)`)\cr
#'   Specifications for MCMC sampling.
#' @param verbose (`logical(1)`)\cr
#'   Whether to print any output to the console or not.
#' @param pool (`logical(1)`)\cr
#'   Whether to fit a pooled BART model for the control data, i.e. not including 
#'   the study/source indicator as a covariate.
#'   For experimental purposes only.
#' @param ... (`Any`)\cr
#'   Additional arguments passed to `dbarts::bart2()`.
#' 
#' @details 
#' The `data` data.frame and the corresponding columns specified via `outcome`, `treatment`, `source` and `covariates` 
#' must meet certain criteria:
#' - `treatment` column must be binary (0, 1), either `numeric()` or `integer()`
#'   - 0 will be interpreted as the control and 1 as the experimental treatment
#' - `source` column must be binary (0, 1), either `numeric()` or `integer()`
#'   - 0 will be interpreted as the current and 1 as the historical trial
#' - all additional categorical covariates must be encoded as factors
#' 
#' @return 
#' `list()`\cr
#' A list with 3 main elements:
#' - `data`: A (quasi-)copy of the data provided to the function.
#' - `bart`: A list containing the two BART models for the treatment and control groups, respectively.
#' - `draws`: A list of draws objects, in particular:
#'   - `y1`: Potential outcomes under the experimental treatment
#'   - `y0`: Potential outcomes under the control treatment
#'   - `cte`: Conditional treatment effects
#'   - `ate`: Average treatment effects
#' 
#' @references 
#' Zhou, Tianjian, and Yuan Ji.
#' „Incorporating External Data into the Analysis of Clinical Trials via Bayesian Additive Regression Trees“.
#' Statistics in Medicine 40, Nr. 28 (2021): 6421–42. https://doi.org/10.1002/sim.9191.
#' 
#' Chipman, Hugh A., Edward I. George, and Robert E. McCulloch.
#' „BART: Bayesian additive regression trees“.
#' The Annals of Applied Statistics 4, Nr. 1 (2010): 266–98. https://doi.org/10.1214/09-AOAS285.
#' 
#' @export
borrow_bart = function(data, outcome, treatment, source, covariates,
                       family = binomial(link = "probit"),
                       num_warmup = 2500L, num_posterior = 2500L, num_chains = 4L,
                       seed = 123,
                       verbose = TRUE,
                       pool = FALSE,
                       ...) {
  # Assertions / sanity checks
  chk$assert_data_frame(data)
  chk$assert_string(outcome)
  chk$assert_subset(outcome, colnames(data))
  chk$assert_string(treatment)
  chk$assert_subset(treatment, colnames(data))
  chk$assert_string(source)
  chk$assert_subset(source, colnames(data))
  chk$assert_character(covariates)
  chk$assert_subset(covariates, colnames(data))
  assert_binary(data[[treatment]])
  chk$assert_class(family, "family")
  chk$assert_subset(family$family, c("gaussian", "binomial"))
  if (family$family == "binomial") {
    # dbarts only implements probit link
    if (family$link == "logit") stop("Must use probit link for binomial models.")
    # Family is inferred in bart function and not determined explicitly; therefore, double-check
    assert_binary(data[[outcome]])
  }
  chk$assert_int(num_warmup, lower = 1)
  chk$assert_int(num_posterior, lower = 1)
  chk$assert_int(num_chains, lower = 1)
  chk$assert_int(seed)
  chk$assert_flag(verbose)
  chk$assert_flag(pool)
  
  
  # Data preparation
  dt = as.data.table(data)
  dt_trt = dt[treatment == 1, , env = list(treatment = treatment)]
  dt_ctrl = dt[treatment == 0, , env = list(treatment = treatment)]
  dt_trial = dt[source == 0, , env = list(source = source)]
  
  
  # Formulas for different BART models
  formula_trt = as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))
  if (pool) {
    formula_ctrl = formula_trt
  } else {
    formula_ctrl = as.formula(paste(outcome, "~", paste(c(covariates, source), collapse = " + ")))
  }
  
  # Treatment BART
  if (verbose) cat("Running MCMC for treatment group...\n\n")
  bart_trt = bart2(
    formula_trt, data = dt_trt,
    # Potential outcomes Y(1)
    test = dt_trial,
    # Match with functions from BART package
    k = 2,
    # MCMC specs
    n.burn = num_warmup, n.samples = num_posterior, n.chains = num_chains,
    # Misc
    verbose = verbose, seed = seed,
    ...
  )
  if (verbose) cat("\n\n")
  
  # Control BART
  if (verbose) cat("Running MCMC for control group...\n\n")
  bart_ctrl = bart2(
    formula_ctrl, data = dt_ctrl,
    # Potential outcomes Y(0)
    test = dt_trial,
    # Match with functions from BART package
    k = 2,
    # MCMC specs
    n.burn = num_warmup, n.samples = num_posterior, n.chains = num_chains,
    # Misc
    verbose = verbose, seed = seed,
    ...
  )
  if (verbose) cat("\n\n")
  
  # For binary data, need to apply response function
  if (family$family == "binomial") {
    bart_trt$yhat.test = pnorm(bart_trt$yhat.test)
    bart_ctrl$yhat.test = pnorm(bart_ctrl$yhat.test)
  }
  
  
  # If 1 chain only: need additional dimension for conversion to draws object
  if (num_chains == 1L) {
    dim(bart_trt$yhat.test) = c(1, dim(bart_trt$yhat.test))
    dim(bart_ctrl$yhat.test) = c(1, dim(bart_ctrl$yhat.test))
  }
  
  # Potential outcome draws
  draws_y1 = as_draws(aperm(bart_trt$yhat.test, c(2, 1, 3)))
  draws_y0 = as_draws(aperm(bart_ctrl$yhat.test, c(2, 1, 3)))
  
  # Conditional treatment effect draws
  draws_cte = draws_y1 - draws_y0
  
  # Average treatment effect draws (obtained using Bayesian bootstrap)
  draws_ate = apply(draws_cte, 1:2, function(x) {
    bb_weights = rdirichlet(n = 1, alpha = rep(1, dim(draws_cte)[3]))
    sum(bb_weights * x)
  })
  # Need to add dimension for conversion to draws object
  dim(draws_ate) = c(dim(draws_ate), 1)
  draws_ate = as_draws(draws_ate)
  
  
  # Final output object
  out = list(
    data = dt,
    bart = list(treatment = bart_trt, control = bart_ctrl),
    draws = list(y1 = draws_y1, y0 = draws_y0, cte = draws_cte, ate = draws_ate)
  )
  class(out) = "borrow_bart"
  
  return(out)
}


# Draw random numbers from dirichlet distribution
rdirichlet = function(n, alpha) {
  x = matrix(rgamma(length(alpha) * n, alpha), ncol = length(alpha))
  x = x / rowSums(x)
  if (n == 1) drop(x) else x
}


