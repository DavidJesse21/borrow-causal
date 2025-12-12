box::use(
  data.table[...],
  stats[glm, binomial, gaussian, fitted, var, rbeta, as.formula, quantile],
  partitions[setparts],
  chk = checkmate
)

box::use(
  ./utils[mapvalues],
  ./utils[assert_binary]
)


#' Bayesian dynamic borrowing using propensity score weighted Multi-source Exchangeability Models
#' 
#' @description
#' This function implements the propensity score weighted multi-source exchangeability model
#' as proposed by Wei et al. (see references below).
#' They combine the multi-source exchangeability modelling framework by Kaizer et al. (see reference below) 
#' with the propensity score methodology for borrowing historical control data.
#' 
#' @details
#' Currently, this implements the PS-MEM method for normal and binary outcome data.
#' For normal data, a flat prior is assumed and there is no option to specify a different kind of prior.
#' For binary data, a beta prior is used with default (hyper-)parameters inducing an uninformative prior.
#' These can in principle be changed, however, in the current implementation they are used for both, 
#' the treatment and control group, so we suggest to leave the parameters at their defaults.
#' Also note that, in the binary case, while the posterior distributions of the two groups can be derived 
#' analytically, the same is not true for the treatment effect.
#' Therefore, Monte Carlo sampling is applied to obtain samples of the risk difference.
#' For the propensity score estimation a logistic regression model is used.
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
#'   Defaults to `NULL` in which case no PS model is estimated and all observations obtain equal weights.
#' @param family (`family`)\cr
#'   A family object (see `?stats::family`).
#'   Possible options are `binomial()` and `gaussian()`.
#' @param prior_delta (`numeric(1)`)\cr
#'   A number used for constructing the prior model probabilities.
#'   A value of 0 corresponds to equal prior probabilities among all candidate models.
#'   Values larger than 0 assign higher prior probabilities to non-exchangeable configurations.
#'   Values smaller than 0 assign higher prior probabilities to exchangeable configurations.
#' @param prior_a,prior_b (`numeric(1)`)\cr
#'   For binomial data, the parameters for the prior beta distribution (used for both control and treatment arm).
#' @param num_mc (`numeric(1)`)\cr
#'   For binomial data, the number of Monte Carlo samples to draw from the posterior distributions of the two treatment groups.
#'   This is required to compute a posterior distribution of the treatment effect.
#'  
#' @return (`list()`)\cr
#'   A `borrow_ps_mem` object, which is a list with three elements:
#'   1. `mem` (`list`): A list containing prior and posterior information about the multi-source exchangeability model.
#'   2. `posterior` (`list()`): A list containing information about the posterior distributions of the control group, treatment group 
#'      as well as the effect of interest.
#'   3. `data` (`data.table()`): The data used for the computation. In addition to the original data,
#'      it has the column `weights` containing the estimated propensity score weights that have been 
#'      used for all further computations.
#'  
#' @references
#' Wei, Wei, Yunxuan Zhang, Satrajit Roychoudhury, and the Alzheimer’s Disease Neuroimaging Initiative.
#' „Propensity Score Weighted Multi-Source Exchangeability Models for Incorporating External Control Data in Randomized Clinical Trials“.
#' Statistics in Medicine 43, Nr. 20 (2024): 3815–29. https://doi.org/10.1002/sim.10158.
#' 
#' Kaizer, Alexander M, Joseph S Koopmeiners, and Brian P Hobbs.
#' „Bayesian hierarchical modeling based on multisource exchangeability“.
#' Biostatistics 19, Nr. 2 (1. April 2018): 169–84. https://doi.org/10.1093/biostatistics/kxx031.
#' 
#' @export
borrow_ps_mem = function(data, outcome, treatment, source, covariates = NULL, family = binomial(),
                         prior_delta = 0, prior_a = 0.5, prior_b = 0.5, num_mc = 5000L) {
  # Assertions / sanity checks
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
  chk$assert_class(family, "family")
  chk$assert_subset(family$family, c("gaussian", "binomial"))
  chk$assert_number(prior_delta)
  chk$assert_number(prior_a, lower = 0)
  chk$assert_number(prior_b, lower = 0)
  chk$assert_number(num_mc, lower = 1L)
  if (any(c("weight", "group") %in% colnames(data))) {
    stop("Columns `weight` and `group` are reserverd names for internal computations ",
         "and must not be contained in `data`.")
  }
  
  
  # data.table for convenience + split by treatment arm
  dt = as.data.table(data)
  dt_trt = dt[treatment == 1, , env = list(treatment = treatment)]
  dt_ctrl = dt[treatment == 0, , env = list(treatment = treatment)]
  
  # Obtain propensity score weights for control arm data
  dt_ctrl = calc_ps_weights(dt_ctrl, source, covariates)
  
  # Initialize MEM object
  mem = init_mem(length(unique(dt[[source]])), prior_delta)
  
  # Calculate posterior model probabilities
  mem$posteriors = calc_posterior_probs(dt_ctrl, outcome, source, "weight", mem, family, prior_a, prior_b)
  
  # Calculate posterior exchangeability probabilities
  mem$exchange_probs = calc_exchange_probs(mem)
  
  # Calculate effective sample sizes
  mem$ess = calc_ess(dt_ctrl, source, "weight", mem$exchange_probs)
  
  # Initialize object for final posterior inference
  posterior = list(
    control = NULL,
    treatment = NULL,
    effect = NULL
  )
  
  # Constant weights needed for internal functions to work
  dt_trt[, weight := 1L]
  
  # Final posterior inference
  if (family$family == "gaussian") {
    # Normal (analytical)
    posterior$control = post_normal(dt_ctrl, outcome, source, "weight", mem)
    posterior$treatment = post_normal(dt_trt, outcome, source, "weight", mem)
    posterior$effect = c(
      mu = posterior$treatment[["mu"]] - posterior$control[["mu"]],
      sigma2 = posterior$treatment[["sigma2"]] + posterior$control[["sigma2"]]
    )
  } else {
    # Binary
    post_ctrl = post_binary(dt_ctrl, outcome, source, "weight", mem)
    post_trt = post_binary(dt_trt, outcome, source, "weight", mem)
    
    # No analytical solution for difference of two Beta RVs -> MC sampling
    x_ctrl = rbeta(num_mc, post_ctrl[["a"]], post_ctrl[["b"]])
    x_trt = rbeta(num_mc, post_trt[["a"]], post_trt[["b"]])
    x_effect = x_trt - x_ctrl
    
    posterior$control = post_ctrl
    posterior$treatment = post_trt
    posterior$effect = list(
      samples = x_effect,
      summary = mc_summary(x_effect)
    )
  }
  
  # Also return data for transparency
  dt_ctrl[, group := NULL]
  dt = rbindlist(list(dt_ctrl, dt_trt))
  
  # Output object
  out = list(
    mem = mem,
    posterior = posterior,
    data = dt
  )
  class(out) = "borrow_ps_mem"
  
  
  return(out)
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
calc_ps_weights = function(dt, source, covariates = NULL) {
  # Special case: estimate no PS model and simply return equal weights
  if (is.null(covariates)) {
    new_dt = copy(dt)
    new_dt[, weight := rep(1L, nrow(new_dt))]
    return(new_dt)
  }
  
  # Make a copy of the data with a "study response", which equals 1 if observation is from current trial and 0 otherwise
  x = copy(dt)
  set(x, j = "s_response", value = fifelse(x[[source]] == 0, 1L, 0L))
  
  # Estimate propensity scores
  ps_formula = as.formula(paste("s_response ~", paste(covariates, collapse = " + ")))
  ps_model = glm(ps_formula, data = x, family = binomial())

  # Calculate propensity score weights
  ps = fitted(ps_model)
  weights = fifelse(x$s_response == 1L, 1, ps / (1 - ps))
  
  # Add propensity score weights to original data
  new_dt = copy(dt)
  new_dt[, weight := weights]
  
  return(new_dt)
}



#' Initialize Multi-source Exchangeability Model (MEM)
#' 
#' @description
#' This function initializes a MEM object by constructing a configuration (exchangeability) matrix 
#' and assigning prior model probabilities.
#' 
#' @param num_studies (`numeric(1)`)\cr
#'   The (positive) number of studies (including the current trial).
#' @param prior_delta (`numeric(1)`)\cr
#'   Used for assigning prior model probabilities.
#'   A value of 0 assigns equal prior probability to all model configurations.
#'   Values larger than 0 assign higher prior probability to non-exchangeable configurations.
#'   Values smaller than 0 assign higher prior probability to exchangeable configurations.
#' 
#' @return (`list()`)\cr
#'   A list with two elements:
#'   * `configs` (`matrix()`):\cr
#'     A matrix denoting the exchangeability configurations.
#'     Rows indicate the configurations, while columns represent the respective studies.
#'   * `priors` (`numeric()`):\cr
#'     A vector containing the prior model probabilities
init_mem = function(num_studies, prior_delta = 0) {
  # Matrix enumerating all exchangeability configurations
  # Columns: studies
  # Rows: exchangeability configuration
  config_mat = as.matrix(t(setparts(num_studies)))
  rownames(config_mat) = paste0("m", 1:nrow(config_mat))
  colnames(config_mat) = paste0("source", 1:ncol(config_mat) - 1)
  
  # Number of configuration and number of parameters within each MEM configuration
  num_configs = nrow(config_mat)
  num_params = apply(config_mat, MARGIN = 1, \(x) length(unique(x)))
  
  # Prior model (exchangeability configuration) probabilities
  model_priors = num_params^prior_delta / sum(num_params^prior_delta)
  
  # Output MEM object
  out = list(configs = config_mat, priors = model_priors)
  
  return(out)
}



#' Calculate the evidence aka marginal likelihood for binary data
#' 
#' @param dt (`data.table()`)\cr
#'   The data. At this time it must contain a weights column.
#' @param outcome (`character(1)`)\cr
#'   The column name of the outcome variable. 
#' @param source (`character(1)`)\cr
#'   The column name of the study source indicator.
#' @param weight (`character(1)`)\cr
#'   The name of the (propensity) weight column.
#' @param a0 (`numeric(1)`)\cr
#'   Alpha parameter for the prior beta distribution.
#' @param b0 (`numeric(1)`)\cr
#'   Beta parameter for the prior beta distribution.
#' @param .log (`logical(1)`)\cr
#'   A flag indicating whether to return the log of the evidence.
#'   
#' @return (`numeric(1)`)\cr
#'   The evidence given the data and the prior.
evidence_binary = function(dt, outcome, source, weight = "weight", a0 = 1, b0 = 1, .log = TRUE) {
  # For binary outcomes the source information is actually not required.
  # Everything boils down to summation.
  n_wt = dt[, sum(weight), env = list(weight = weight)]
  n_events_wt = dt[, sum(weight * outcome), env = list(outcome = outcome, weight = weight)]
  
  a1 = a0 + n_events_wt
  b1 = b0 + n_wt + n_events_wt
  
  out = lbeta(a1, b1) - lbeta(a0, b0)
  
  if (!.log) {
    return(exp(out))
  } else {
    return(out)
  }
}


#' Calculate the evidence aka marginal likelihood for normally distributed data
#' 
#' @param dt (`data.table()`)\cr
#'   The data. At this time it must contain a weight column.
#' @param outcome (`character(1)`)\cr
#'   The column name of the outcome variable. 
#' @param source (`character(1)`)\cr
#'   The column name of the study source indicator.
#' @param weight (`character(1)`)\cr
#'   The name of the (propensity) weight column.
#' @param .log (`logical(1)`)\cr
#'   A flag indicating whether to return the log of the evidence.
#'   
#' @return (`numeric(1)`)\cr
#'   The evidence given the data and the prior.
evidence_normal = function(dt, outcome, source, weight = "weight", .log = TRUE) {
  li_dt = split(dt, by = source)
  
  # Calculate "u" and "v" terms (see Table 1, p.3819)
  x = vapply(li_dt, function(dts) {
    wvar = dts[, weighted_var(outcome, weight), env = list(outcome = outcome, weight = weight)]
    upart = dts[, sum(weight * outcome) / wvar, env = list(outcome = outcome, weight = weight)]
    vpart = dts[, 0.5 * sum(weight) / wvar, env = list(outcome = outcome, weight = weight)]
    return(c(u = upart, v = vpart))
  }, numeric(2)) |>
    rowSums()
  
  t1 = x[["u"]]^2 / (4 * x[["v"]])
  t2 = 0.5 * (log(pi) - log(x[["v"]]))
  out = t1 + t2
  
  if (!.log) {
    return(exp(out))
  } else {
    return(out)
  }
}




#' Calculate posterior model probabilities
#' 
#' @param dt (`data.table()`)\cr
#'   The data.
#' @param outcome (`character(1)`)\cr
#'   The column name of the outcome variable.
#' @param source (`character(1)`)\cr
#'   The column name of the study source indicator.
#' @param weight (`character(1)`)\cr
#'   The name of the (propensity) weight column.
#' @param mem (`list()`)\cr
#'   A `mem` object as returned by `init_mem()`.
#' @param family (`family`)\cr
#'   A family object (see `?stats::family`).
#'   Possible options are `binomial()` and `gaussian()`.
#' @param a0 (`numeric(1)`)\cr
#'   Alpha parameter for the prior beta distribution.
#' @param b0 (`numeric(1)`)\cr
#'   Beta parameter for the prior beta distribution.
#'   
#' @return (`numeric()`)\cr
#'   A (named) vector containing the posterior model probabilities.
calc_posterior_probs = function(dt, outcome, source, weight = "weight", mem, family, a0, b0) {
  # Auxiliary preparations for calculation of log evidences
  configs = mem$configs
  colnames(configs) = sub("^source", "", colnames(configs))
  dt[, group := NA_integer_]
  
  # (Sum of) log evidences (marginal likelihoods) for each model configuration
  log_evidences = apply(configs, MARGIN = 1, function(config) {
    # Assign group/cluster label given the current model/exchangeability pattern
    set(
      dt, j = "group",
      value = mapvalues(dt[[source]], from = as.integer(names(config)), to = unname(config))
    )
    # Split data by this group
    li_dt = split(dt, by = "group")
    # Actual calculations
    if (family$family == "gaussian") {
      evidence_parts = vapply(li_dt, \(x) evidence_normal(x, outcome, source, weight), numeric(1))
      sum(evidence_parts)
    } else {
      evidence_parts = vapply(li_dt, \(x) evidence_binary(x, outcome, soure, weight, a0, b0), numeric(1))
      sum(evidence_parts)
    }
  })
  
  # Sum with log prior model probabilities
  log_prior_evidences = log(mem$priors) + log_evidences
  
  # Final model probabilities
  posterior_probs = exp(log_prior_evidences) / sum(exp(log_prior_evidences))
  
  
  return(posterior_probs)
}



#' Calculate (posterior) exchangeability probabilities
#' 
#' @param mem (`list()`)\cr
#'   A `mem` object already containing the posterior model probabilities.
#' @return (`numeric()`)\cr
#'   A (named) vector containing the posterior exchangeability probabilities for the different 
#'   (historical) studies. 
calc_exchange_probs = function(mem) {
  configs = mem$configs
  probs = mem$posteriors
  idx = apply(configs, MARGIN = 2, \(x) as.integer(x == configs[, 1]))
  exchange_probs = colSums(probs * idx)
  
  return(exchange_probs)
}



#' Obtain posterior distribution for a normal likelihood with flat prior
#' 
#' @param dt (`data.table()`)\cr
#'   The data.
#' @param outcome (`character(1)`)\cr
#'   The column name of the outcome variable.
#' @param source (`character(1)`)\cr
#'   The column name of the study source indicator.
#' @param weight (`character(1)`)\cr
#'   The name of the (propensity) weight column.
#' @param mem (`list()`)\cr
#'   A `mem` object already containing the posterior exchangeability probabilities.
#'  
#' @details
#' For this derivation of the posterior distribution, it must be assumed that the variance of the 
#' likelihood model is known.
#' In fact, it must be estimated, thus we use a small hack by assuming that the estimated variance 
#' equals the unknown true variance.
#'   
#' @return (`numeric(2)`)\cr
#'   A named numeric vector containing the posterior mean (`"mu"`) and the posterior variance (`"sigma2"`).
post_normal = function(dt, outcome, source, weight = "weight", mem) {
  exchange_probs = mem$exchange_probs
  
  vals = dt[, .(n_wt = sum(weight),
                wmean = weighted_mean(outcome, weight),
                wvar = weighted_var(outcome, weight)),
            by = source,
            env = list(weight = weight, outcome = outcome, source = source)]
  setorderv(vals, source)
  vals = as.matrix(vals)
  
  # Posterior precision
  prec = sum((vals[, "n_wt"] / vals[, "wvar"]) * exchange_probs)
  # Posterior mean
  mu = sum(vals[, "n_wt"] * vals[, "wmean"] * exchange_probs / vals[, "wvar"]) / prec
  
  return(c(mu = mu, sigma2 = 1/prec))
}



#' Obtain posterior distribution for a binomial likelihood with beta prior
#' 
#' @param dt (`data.table()`)\cr
#'   The data.
#' @param outcome (`character(1)`)\cr
#'   The column name of the outcome variable.
#' @param source (`character(1)`)\cr
#'   The column name of the study source indicator.
#' @param weight (`character(1)`)\cr
#'   The name of the (propensity) weight column.
#' @param mem (`list()`)\cr
#'   A `mem` object already containing the posterior exchangeability probabilities.
#' @param a0 (`numeric(1)`)\cr
#'   Alpha parameter for the prior beta distribution.
#' @param b0 (`numeric(1)`)\cr
#'   Beta parameter for the prior beta distribution.
#'   
#' @return (`numeric(2)`)\cr
#'   A named numeric vector containing the parameters of the beta posterior distribution.
post_binary = function(dt, outcome, source, weight = "weight", mem, a0 = 1, b0 = 1) {
  exchange_probs = mem$exchange_probs

  vals = dt[, .(n_wt = sum(weight),
                n_events_wt = sum(outcome * weight)),
            by = source,
            env = list(weight = weight, outcome = outcome, source = source)]
  setorderv(vals, source)
  vals = as.matrix(vals)
  
  a1 = a0 + sum(vals[, "n_events_wt"] * exchange_probs)
  b1 = b0 + sum((vals[, "n_wt"] - vals[, "n_events_wt"]) * exchange_probs)
  
  return(c(a = a1, b = b1))
}



#' Calculate the effective sample size of the (external) control data
#' 
#' @param dt (`data.table()`)\cr
#'   A data.table containing the control arm data as well as their associated propensity 
#'   score weight.
#' @param source (`character(1)`)\cr
#'   A string specifying the column name for the source/study indicator variable.
#' @param weight (`character(1)`)\cr
#'   The name of the (propensity) weight column.
#' @param exchange_probs (`numeric()`)\cr
#'   A vector of exchangeability probabilites.
#' 
#' @return (`numeric()`)\cr
#' A vector containing the effective sample sizes of the different control arm data sources.  
#'   
#' @note
#' For the current/internal control arm the "effective" sample size is also calculated 
#' for consistency reasons.
calc_ess = function(dt, source, weight = "weight", exchange_probs) {
  idx = dt[, sort(unique(source)), env = list(source = source)]
  
  ess = vapply(idx, function(i) {
    # Internal control arm, ESS is simply the sample size
    if (i == 0) {
      dt[source == i, .N, env = list(source = source)]
    } else {
      # Propensity score weighted external controls
      dt[source == i,
         sum(weight * exchange_probs[[paste0("source", i)]]),
         env = list(source = source, weight = weight)]
    }
  }, numeric(1))
  
  names(ess) = paste0("source", idx)
  
  return(ess)
}


#' Small but important helpers
weighted_mean = function(x, weight) {
  sum(weight * x) / sum(weight)
}

weighted_var = function(x, weight) {
  mu = weighted_mean(x, weight)
  sum(weight * (x - mu)^2) / (sum(weight) - 1)
}

mc_summary = function(x, probs = c(0.025, 0.5, 0.975)) {
  x1 = mean(x)
  x2 = quantile(x, probs = probs)
  c(Mean = x1, x2)
}


