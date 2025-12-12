box::use(
  R6[R6Class],
  mvtnorm[rmvnorm],
  stats[rnorm],
  data.table[...],
  chk = checkmate
)


#' Simulator object
#' 
#' @description
#' This simulator object contains all parameter values required for the simulation model
#' to be fully specified.
#' They can be set/changed when initializing a new instance of the object or via the 
#' `set_params()` method.
#' A single data set can be generated via the `gen_data()` method.
#' A number of replicates can be generated via the `simulate()` method.
#' 
#' @export
Simulator = R6Class(
  "Simulator",
  
  public = list(
    
    #' @field params (`list()`)\cr
    #'   List of all parameters of the simulation model.
    params = list(
      # Sample sizes
      n_trt = 60,
      n_ctrl = 30,
      n_hist = 30,
      # Covariate moments
      x_mean = rep(8, 3),
      x_cov = matrix(c(6, 2, 0, 2, 6, 0, 0, 0, 6), nrow = 3, byrow = TRUE),
      # Regression coefficients
      b_0 = 2,
      b_x = c(1, 1, 0),
      b_trt = 0,
      # Observed and unobserved variable shifts
      shift_x2 = 0,
      shift_x3 = 0,
      shift_u = 0,
      # Nuisance parameters
      sd_u = 4,
      sd_err = 7
    ),
    
    
    #' @description
    #' Create a new simulator object.
    #' 
    #' @param ...,li_params See `set_params()` method.
    initialize = function(..., li_params = NULL) {
      self$set_params(..., li_params = li_params)
      return(invisible(self))
    },
    
    
    #' @description
    #' Set the parameters of the simulation model.
    #' @param ... (`Any`)\cr
    #'   Named arguments relating to the model parameters to be set.
    #' @param li_params (`list()`)\cr
    #'   Provide a named list of parameters alternatively to using `...`.
    #'
    #' @details
    #' The following parameters can/must be set:
    #' - `n_trt`, `n_ctrl`, `n_hist` (`numeric(1)`):
    #'   The sample size for the internal treatment, internal control, and external control group, respectively.
    #' - `x_mean` (`numeric(3)`): Mean vector for joint normal distribution of covariates.
    #' - `x_cov` (`matrix(numeric(9), nrow = 3, ncol = 3)`):
    #'    Covariance matrix for joint normal distribution of covariates.
    #' - `b_0` (`numeric(1)`): Intercept of linear regression model.
    #' - `b_x` (`numeric(3)`): Regression coefficients for the three covariates.
    #' - `b_trt` (`numeric(1)`): Regression coefficient for the treatment effect.
    #' - `shift_x2`, `shift_x3`, `shift_u` (`numeric(1)`):
    #'   Mean shift in the external control population relative to the current trial population w.r.t.
    #'   the covariates `x2` and `x3`, and the latent variable `u`.
    #' - `sd_u`, `sd_err` (`numeric(1)`): Standard deviation of latent variable `u` and 
    #'   the random error term.
    set_params = function(..., li_params = NULL) {
      # Make sure only ... or li_params is used
      if (...length() > 0 && !is.null(li_params)) {
        stop("Either use `...` or `li_params` but not both.")
      } else if (...length() == 0 && is.null(li_params)) {
        return(invisible(self))
      } else if (...length() > 0) {
        li_params = list(...)
      } else {
        li_params = li_params
      }
      
      # Make sure supplied input is valid
      chk$assert_subset(names(li_params), names(self$params), empty.ok = FALSE)
      for (param in names(li_params)) {
        private$assertion_list[[param]](li_params[[param]])
      }
      # Set new values
      self$params[names(li_params)] = li_params
      
      return(invisible(self))
    },
    
    
    #' @description
    #' Simulate data sets.
    #' @param num_sims (`numeric(1)`)\cr
    #'   The number of data sets to be simulated.
    simulate = function(num_sims) {
      for (param in names(self$params)) {
        if (any(is.na(self$params[[param]]))) {
          stop(sprintf("Parameter `%s` needs to be set.", param))
        }
      }
      chk$assert_count(num_sims)
      
      lapply(1:num_sims, \(i) self$gen_data(check_params = FALSE))
    },
    
    
    #' @description
    #' Generate a single data set.
    #' @param check_params (`logical(1)`)\cr
    #'   Whether to (explicitly) check if all parameters are set before trying to
    #'   generate a data set.
    gen_data = function(check_params = TRUE) {
      if (check_params) {
        for (param in names(self$params)) {
          if (any(is.na(self$params[[param]]))) {
            stop(sprintf("Parameter `%s` needs to be set.", param))
          }
        }
      }
      
      dt_trt = private$.gen_data(
        n = self$params$n_trt,
        trt = 1L,
        x_mean = self$params$x_mean,
        u_mean = 0,
        source = 0L
      )
      
      dt_ctrl = private$.gen_data(
        n = self$params$n_ctrl,
        trt = 0L,
        x_mean = self$params$x_mean,
        u_mean = 0,
        source = 0L
      )
      
      dt_hist = private$.gen_data(
        n = self$params$n_hist,
        trt = 0L,
        x_mean = self$params$x_mean + c(0, self$params$shift_x2, self$params$shift_x3),
        u_mean = self$params$shift_u,
        source = 1L
      )
      
      dt = rbindlist(list(dt_trt, dt_ctrl, dt_hist))
      return(dt)
    }
  ),
  
  
  private = list(
    # Internal data generating function
    .gen_data = function(n, trt, x_mean, u_mean, source) {
      x = rmvnorm(n, x_mean, self$params$x_cov)
      linpred = self$params$b_0 + drop(cbind(x, rep(trt, n)) %*% c(self$params$b_x, self$params$b_trt))
      u = rnorm(n, u_mean, self$params$sd_u)
      err = rnorm(n, 0, self$params$sd_err)
      y = linpred + u + err
      
      dt = cbind(x, rep(trt, n), rep(source, n), y)
      colnames(dt) = c("x1", "x2", "x3", "trt", "source", "y")
      dt = as.data.table(dt)
      dt[, `:=`(trt = as.integer(trt),
                source = as.integer(source))]
      
      return(dt[])
    },
    
    # List with assertion functions for each parameter
    assertion_list = list(
      # Sample sizes
      n_trt = \(n_trt) chk$assert_count(n_trt),
      n_ctrl = \(n_ctrl) chk$assert_count(n_ctrl),
      n_hist = \(n_hist) chk$assert_count(n_hist),
      # Covariate moments
      x_mean = \(x_mean) chk$assert_numeric(x_mean, len = 3),
      x_cov = \(x_cov) chk$assert_matrix(x_cov, mode = "numeric", nrows = 3, ncols = 3),
      # Regression coefficients
      b_0 = \(b_0) chk$assert_number(b_0),
      b_x = \(b_x) chk$assert_numeric(b_x, len = 3),
      b_trt = \(b_trt) chk$assert_number(b_trt),
      # Observed and unobserved variable shifts
      shift_x2 = \(shift_x2) chk$assert_number(shift_x2),
      shift_x3 = \(shift_x3) chk$assert_number(shift_x3),
      shift_u = \(shift_u) chk$assert_number(shift_u),
      # Nuisance parameters
      sd_u = \(sd_u) chk$assert_number(sd_u, lower = 0),
      sd_err = \(sd_err) chk$assert_number(sd_err, lower = 0)
    )
  )
)


#' Simulation scenarios
#' 
#' @param drop_constants (`logical(1)`)\cr
#'   Whether to drop columns with a single value only.
#' 
#' @export
make_scenarios = function(drop_constants = TRUE) {
  scenarios = CJ(
    # Sample sizes
    n_trt = 60,
    n_ctrl = 30,
    n_hist = c(1, 2, 4) * 30,
    # Treatment effect
    b_trt = c(0, 3.5, 4.5),
    # x2 prognostic for Y -> observed confounding
    shift_x2 = seq(-3, 3, by = 1),
    # x3 not prognostic for Y -> no confounding, i.e. irrelevant/misleading information
    # (Therefore, no need to cover both directions)
    shift_x3 = seq(0, 3, by = 1),
    # u = unobserved/latent confounding
    shift_u = seq(-3, 3, by = 1)
  )
  
  # Only consider one direction of confounding (combination of observed and unobserved)
  scenarios = scenarios[sign(shift_x2 * shift_u) >= 0]
  
  # Drop columns with constant values if requested
  if (drop_constants) {
    for (j in c("n_trt", "n_ctrl")) set(scenarios, j = j, value = NULL)
  }
  
  # Add identifier
  scenarios[, scenario.id := 1:.N]
  setcolorder(scenarios, neworder = "scenario.id")
  setkey(scenarios, scenario.id)
  
  return(scenarios[])
}


# Test simulation model ----

if (is.null(box::name())) {
  box::use(
    testthat[...],
    data.table[...],
    stats[var, cov]
  )
  
  
  test_that("Data structur is correct", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 110, n_ctrl = 100, n_hist = 90)
    dt = simulator$gen_data()
    expect_equal(colnames(dt), c("x1", "x2", "x3", "trt", "source", "y"))
    expect_equal(nrow(dt), 300)
    expect_equal(dt[, .N, by = .(trt, source)]$N, c(110, 100, 90))
  })
  
  
  test_that("Distributions of x1, x2, and x3 are correct", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 10, n_ctrl = 1000, n_hist = 10)
    dt = simulator$gen_data()
    dt = dt[trt == 0 & source == 0]
    
    expect_equal(dt[, mean(x1)], simulator$params$x_mean[1], tolerance = 1)
    expect_equal(dt[, var(x1)], diag(simulator$params$x_cov)[1], tolerance = 1)
    
    expect_equal(dt[, mean(x2)], simulator$params$x_mean[2], tolerance = 1)
    expect_equal(dt[, var(x2)], diag(simulator$params$x_cov)[2], tolerance = 1)
    
    expect_equal(dt[, mean(x3)], simulator$params$x_mean[3], tolerance = 1)
    expect_equal(dt[, var(x3)], diag(simulator$params$x_cov)[3], tolerance = 1)
    
    expect_equal(dt[, cov(x1, x2)], simulator$params$x_cov[1, 2], tolerance = 0.5)
    expect_equal(dt[, cov(x1, x3)], simulator$params$x_cov[1, 3], tolerance = 0.5)
    expect_equal(dt[, cov(x2, x3)], simulator$params$x_cov[2, 3], tolerance = 0.5)
  })
  
  
  test_that("Marginal distribution of y is correct", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 10, n_ctrl = 2000, n_hist = 10)
    dt = simulator$gen_data()
    dt = dt[trt == 0]
    expect_equal(dt[, mean(y)], 18, tolerance = 0.5)
    expect_equal(dt[, var(y)], 81, tolerance = 2)
  })
  
  
  test_that("Treatment effect gets added", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 1000, n_ctrl = 1000, n_hist = 10, b_trt = 4)
    dt = simulator$gen_data()
    dt = dt[source == 0]
    expect_equal(dt[, mean(y), by = trt][, diff(rev(V1))], 4, tolerance = 0.5)
  })
  
  
  test_that("Covariate shifts work", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 2000, n_ctrl = 2000, shift_x2 = -1.5, shift_x3 = -1.5)
    dt = simulator$gen_data()
    dt = dt[trt == 0]
    
    # Test the shifts themselves
    expect_equal(dt[, mean(x1), by = source][, diff(rev(V1))], 0, tolerance = 0.25)
    expect_equal(dt[, mean(x2), by = source][, diff(rev(V1))], 1.5, tolerance = 0.25)
    expect_equal(dt[, mean(x3), by = source][, diff(rev(V1))], 1.5, tolerance = 0.25)
    
    # Test if the shifts are correctly propagated to the outcome
    # Note: only shift of x2 impacts y
    expect_equal(dt[, mean(y), by = source][, diff(rev(V1))], 1.5, tolerance = 0.5)
  })
  
  
  test_that("Unobserved confounding works", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 10, n_ctrl = 1000, n_hist = 1000, shift_u = -3)
    dt = simulator$gen_data()
    dt = dt[trt == 0]
    
    expect_equal(dt[, mean(y), by = source][, diff(rev(V1))], 3, tolerance = 0.5)
  })
  
  
  test_that("Simulating multiple data sets works", {
    set.seed(123)
    simulator = Simulator$new(n_trt = 60, n_ctrl = 30, n_hist = 30)
    li_dt = simulator$simulate(3)
    expect_equal(length(li_dt), 3)
    expect_true(all(vapply(li_dt, nrow, numeric(1)) == 120))
  })
  
  
  test_file(box::file())
}

