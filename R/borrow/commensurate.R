box::use(
  data.table[...],
  stats[glm, binomial, gaussian, predict, as.formula, fitted],
  psb = psborrow2,
  chk = checkmate
)

box::use(
  ./utils[assert_binary]
)


#' Initialize model using a propensity-score weighted commensurate prior
#' 
#' @description
#' This function initializes a model for employing the PS-weighted commensurate prior.
#' In the end, an `Analysis` object from the `psborrow2` package will be created and, by that, a Stan model will be compiled.
#' No MCMC is run at this point.
#' 
#' 
#' @param data_init (`data.frame()`)\cr
#'   An (initial) data.frame containing the combined current and historical trial data.
#' @param outcome (`character(1)`)\cr
#'   A string specifying the column name for the outcome variable.
#' @param treatment (`character(1)`)\cr
#'   A string specifying the column name for the treatment variable.
#' @param source (`character(1)`)\cr
#'   A string specifying the column name for the source/study indicator variable.
#' @param covariates (`character()`)\cr
#'   A character vector specifying the column names of covariates to be used in the models.
#'   Defaults to `NULL` in which case no PS model is estimated and all observations obtain equal weights.
#' @param family (`family`)\cr
#'   A family object (see `?stats::family`).
#'   Possible options are `binomial()` and `gaussian()`.
#' @param prior_tau (`S4`)\cr
#'   Prior for the commensurability parameter.
#'   Needs to be created with a suitable function from the `psborrow2` package.
#' @param prior_baseline (`S4`)\cr
#'   Prior for the intercept of the external control group.
#'   Needs to be created with a suitable function from the `psborrow2` package.
#' @param prior_std_dev (`S4`)\cr
#'   Prior for the residual standard deviation in case of Gaussian data.
#'   Needs to be created with a suitable function from the `psborrow2` package.
#' @param prior_trt (`S4`)\cr
#'   Prior for the treatment effect.
#'   Needs to be created with a suitable function from the `psborrow2` package.
#'  
#' @return (`Analysis`)\cr
#'   An `Analysis` object from the `psborrow2` package.
#'  
#' @references
#' Wang, Xi, Leah Suttner, Thomas Jemielita, und Xiaoyun Li.
#' „Propensity score-integrated Bayesian prior approaches for augmented control designs: a simulation study“.
#' Journal of Biopharmaceutical Statistics 32, Nr. 1 (2. Januar 2022): 170–90. https://doi.org/10.1080/10543406.2021.2011743.
#' 
#' Hobbs, Brian P., Bradley P. Carlin, Sumithra J. Mandrekar, und Daniel J. Sargent.
#' „Hierarchical Commensurate and Power Prior Models for Adaptive Incorporation of Historical Information in Clinical Trials“.
#' Biometrics 67, Nr. 3 (1. September 2011): 1047–56. https://doi.org/10.1111/j.1541-0420.2011.01564.x.
#' 
#' @export
init_commensurate = function(data_init, outcome, treatment, source, covariates = NULL,
                             family = gaussian(),
                             prior_tau = psb$prior_gamma(alpha = 0.001, beta = 0.001),
                             prior_baseline = psb$prior_normal(0, 100),
                             prior_std_dev = psb$prior_half_cauchy(0, 5),
                             prior_trt = psb$prior_normal(0, 100)) {
  # Assertions / sanity checks
  chk$assert_data_frame(data_init)
  chk$assert_string(outcome)
  chk$assert_subset(outcome, colnames(data_init))
  chk$assert_string(treatment)
  chk$assert_subset(treatment, colnames(data_init))
  chk$assert_string(source)
  chk$assert_subset(source, colnames(data_init))
  chk$assert_character(covariates, null.ok = TRUE)
  chk$assert_subset(covariates, colnames(data_init))
  assert_binary(data_init[[treatment]])
  chk$assert_class(family, "family")
  chk$assert_subset(family$family, c("gaussian", "binomial"))
  assert_prior(prior_tau)
  assert_prior(prior_baseline)
  assert_prior(prior_std_dev)
  assert_prior(prior_trt)
  
  # Add propensity score weights to data
  data_init = calc_ps_weights(data_init, treatment, source, covariates)
  
  # Create initial data matrix for creation of Stan model object
  bdb_data_init = psb$create_data_matrix(
    data_init,
    outcome = outcome,
    trt_flag_col = treatment,
    ext_flag_col = source,
    weight_var = "weight"
  )
  
  # Commensurate prior specification
  bdb_prior = psb$borrowing_hierarchical_commensurate(
    ext_flag_col = source,
    tau_prior = prior_tau
  )
  
  # Outcome specification
  if (family$family == "gaussian") {
    bdb_outcome = psb$outcome_cont_normal(
      continuous_var = outcome,
      # Prior for intercept of external control arm
      baseline_prior = prior_baseline,
      # Residual standard deviation prior
      std_dev_prior = prior_std_dev,
      weight_var = "weight"
    )
  } else if (family$family == "binomial") {
    bdb_outcome = psb$outcome_bin_logistic(
      binary_var = outcome,
      baseline_prior = prior_baseline,
      weight_var = "weight"
    )
  } else {
    NULL
  }
  
  # Treatment variable/effect specifications
  bdb_trt = psb$treatment_details(
    trt_flag_col = treatment,
    trt_prior = prior_trt
  )
  
  # Compile Stan model
  bdb_ana = psb$create_analysis_obj(
    data_matrix = bdb_data_init,
    outcome = bdb_outcome,
    borrowing = bdb_prior,
    treatment = bdb_trt
  )
  
  return(bdb_ana)
}


#' Run MCMC for commensurate prior model
#' 
#' @param model (`S4`)\cr
#'   An `Analysis` object created with the `psborrow2` package.
#' @param data (`data.frame()`)\cr
#'   An (initial) data.frame containing the combined current and historical trial data.
#' @param outcome (`character(1)`)\cr
#'   A string specifying the column name for the outcome variable.
#' @param treatment (`character(1)`)\cr
#'   A string specifying the column name for the treatment variable.
#' @param source (`character(1)`)\cr
#'   A string specifying the column name for the source/study indicator variable.
#' @param covariates (`character()`)\cr
#'   A character vector specifying the column names of covariates to be used in the models.
#'   Defaults to `NULL` in which case no PS model is estimated and all observations obtain equal weights.
#' @param num_warmup,num_posterior,num_chains,seed (`numeric(1)`)\cr
#'   Arguments for MCMC sampling passed to `psborrow2::mcmc_sample()`.
#' @param verbose (`logical(1)`)\cr
#'   Passed to `psborrow2::mcmc_sample()` to control console output during MCMC sampling.
#' @param ... (`Any`)\cr
#'   Additional arguments passed to `psborrow2::mcmc_sample()` or to `$sample()` method of `cmdstanr` model, respectively.
#'   
#' @return (`CmdStanMCMC`)\cr
#'   A `CmdStanMCMC` object containing the posterior draws.
#'   All typical methods and functions from the Stan ecosystem in R are applicable.
#'   
#' @export
borrow_commensurate = function(model, data,
                               outcome, treatment, source, covariates = NULL,
                               num_warmup = 2500, num_posterior = 2500, num_chains = 4L,
                               seed = 123,
                               verbose = FALSE,
                               ...) {
  # Assertions / sanity checks
  assert_analysis(model)
  chk$assert_data_frame(data)
  chk$assert_string(outcome)
  chk$assert_subset(outcome, colnames(data))
  chk$assert_string(treatment)
  chk$assert_subset(treatment, colnames(data))
  chk$assert_string(source)
  chk$assert_subset(source, colnames(data))
  chk$assert_character(covariates, null.ok = TRUE)
  chk$assert_subset(covariates, colnames(data))
  assert_binary(data[[treatment]])
  chk$assert_int(num_warmup, lower = 1)
  chk$assert_int(num_posterior, lower = 1)
  chk$assert_int(num_chains, lower = 1)
  chk$assert_int(seed)
  chk$assert_flag(verbose)
  
  # Create new data to which the model should be fitted (including propensity score weights)
  new_data = calc_ps_weights(data, treatment, source, covariates)
  bdb_data = psb$create_data_matrix(
    new_data,
    outcome = outcome,
    trt_flag_col = treatment,
    ext_flag_col = source,
    weight_var = "weight"
  )
  assert_bdb_data(bdb_data, old_data = model@data_matrix)
  model@data_matrix = bdb_data
  
  # MCMC sampling
  mcmc = psb$mcmc_sample(
    model,
    iter_warmup = num_warmup,
    iter_sampling = num_posterior,
    chains = num_chains,
    verbose = verbose,
    ...
  )
  
  return(mcmc)
}



#' Calculate propensity score weights
#' 
#' @description
#' Fits a generalized linear model with logit link function to calculate study propensity scores 
#' and corresponding weights.
#' 
#' @param dt (`data.table()`)\cr
#'   The data used for propensity score estimation.
#' @param source (`character(1)`)\cr
#'   The column name of the study source indicator.
#' @param covariates (`character()`)\cr
#'   A list of variable names contained in `dt` to be included in the propensity score model.
#'   
#' @return (`data.table()`)\cr
#'   A data.table with the column `weight` attached containing the propensity score weights.
#'   Note that the row ordering of the might be different than that of the inputted data.
calc_ps_weights = function(dt, treatment, source, covariates = NULL) {
  # Special case: estimate no PS model and simply return equal weights
  if (is.null(covariates)) {
    new_dt = copy(dt)
    new_dt[, weight := 1L]
    return(new_dt)
  }
  
  # Otherwise: first split data by treatment (active treatment arm not involved in PS estimation)
  li_dt = split(dt, by = treatment)
  
  dt0 = li_dt[["0"]]
  dt1 = li_dt[["1"]]
  
  # All active treatment subjects: weight of 1
  dt1[, weight := 1L]
  
  # Estimate PS model
  dt0_ps = copy(dt0)
  set(dt0_ps, j = "s_response", value = fifelse(dt0_ps[[source]] == 0, 1L, 0L))
  ps_formula = as.formula(paste("s_response ~", paste(covariates, collapse = " + ")))
  ps_model = glm(ps_formula, data = dt0_ps, family = binomial())
  
  # Calculate PS weights
  psw = fitted(ps_model)
  psw = psw / (1 - psw)
  dt0[, weight := fifelse(source == 0, 1, psw), env = list(source = source)]
  
  out = rbindlist(list(dt0, dt1))
  return(out)
}



# Custom assertion functions ----

# Assertion function for prior objects from `psborrow2` package
checkPrior = function(x) {
  check1 = startsWith(class(x), "Prior")
  check2 = (attr(class(x), "package") == "psborrow2")
  res = check1 && check2
  
  if (isFALSE(res)) {
    return("Must be a prior from the `psborrow2` package")
  } else {
    return(TRUE)
  }
}
assert_prior = chk$makeAssertionFunction(checkPrior)


# Assertion function for analysis objects from `psborrow2` package
checkAnalysis = function(x) {
  check1 = inherits(x, "Analysis")
  check2 = (attr(class(x), "package") == "psborrow2")
  res = check1 && check2
  
  if (isFALSE(res)) {
    return("Must be an `Analysis` object from the `psborrow2` package")
  } else {
    return(TRUE)
  }
}
assert_analysis = chk$makeAssertionFunction(checkAnalysis)


# Double check that new data is valid
checkBdbData = function(new_data, old_data) {
  check1 = is.matrix(new_data)
  check2 = (ncol(new_data) == ncol(old_data))
  check3 = all.equal(colnames(new_data), colnames(old_data))
  res = all(check1, check2, check3)
  
  if (isFALSE(res)) {
    return("Data does not match the model specification")
  } else {
    return(TRUE)
  }
}
assert_bdb_data = chk$makeAssertionFunction(checkBdbData)
