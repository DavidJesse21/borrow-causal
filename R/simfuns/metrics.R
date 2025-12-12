box::use(
  stats[var, na.omit],
  data.table[between]
)


#' Performance metrics
#' 
#' @param x (`numeric()`)\cr
#'   The point estimates of the treatment effect.
#' @param true_val (`numeric(1)`)\cr
#'   True value of the treatment effect.
#' @param na.rm (`logical(1)`)\cr
#'   Whether to remove `NA`s from the calculations.
#' @param ci_lower,ci_upper (`numeric()`)\cr
#'   Lower and upper bounds of the confidence/credible intervals.
#' @param pval (`numeric(1)`)\cr
#'   P-values from the simulations.
#' @param level (`numeric(1)`)\cr
#'   Nominal level for (one-sided) hypothesis testing.
#' @param prop (`logical(1)`)\cr
#'   For calculating the amount/number of non-existing results:
#'   Whether to calculate the relative (percentage) amount (`TRUE`) or 
#'   the absolute number (`FALSE`).
#' 
#' @returns (`numeric(1)`)\cr
#'   The respective performance metric.
#' 
#' @name metrics


#' @rdname metrics
#' @export
calc_bias = function(x, true_val, na.rm = TRUE) {
  mean(x - true_val, na.rm = na.rm) 
}

#' @rdname metrics
#' @export
calc_emp_se = function(x, na.rm = TRUE) {
  sqrt(var(x, na.rm = na.rm))
}

#' @rdname metrics
#' @export
calc_mse = function(x, true_val, na.rm = TRUE) {
  if (na.rm) {
    if (length(true_val) > 1) true_val = true_val[!is.na(x)]
    x = na.omit(x)
  }
  
  sum((x - true_val)^2) / length(x)
}

#' @rdname metrics
#' @export
calc_ci_coverage = function(ci_lower, ci_upper, true_val, na.rm = TRUE) {
  mean(between(true_val, ci_lower, ci_upper, incbounds = TRUE), na.rm = na.rm)
}

#' @rdname metrics
#' @export
calc_ci_mean_width = function(ci_lower, ci_upper, na.rm = TRUE) {
  mean(ci_upper - ci_lower, na.rm = na.rm)
}

#' @rdname metrics
#' @export
calc_rejection_rate = function(pval, level = 0.025, na.rm = TRUE) {
  mean(pval <= level, na.rm = na.rm)
}

#' @rdname metrics
#' @export
calc_num_na = function(x, prop = TRUE) {
  num_na = sum(is.na(x))
  if (prop) {
    return(num_na / length(x))
  } else {
    return(num_na)
  }
}
