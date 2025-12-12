options(box.path = "R")

box::use(
  simfuns/gen_data[Simulator, make_scenarios],
  borrow/commensurate[init_commensurate],
  simfuns/algos[wrapper_all_algos]
)

box::use(
  data.table[...],
  withr[with_seed],
  chk = checkmate,
  stats[gaussian],
  parallel[mclapply],
  parallelly[availableCores],
  fs
)


#' Run simulations for one scenario
#' 
#' @param row_params (`data.table()`)\cr
#'   A one-row data.table containing the parameters for the simulation model/scenario.
#' @param num_sims (`numeric(1)`)\cr
#'   The number of simulation repetitions.
#' @param constants (`list()`)\cr
#'   A named list of constant values (e.g. significance level or algorithm specifications) 
#'   passed to `wrapper_all_algos()`.
#' @param num_cores (`numeric(1)`)\cr
#'   Number of cores across which to parallelize the computations on the generated data sets.
#' @param stan_output_dir (`character(1)`)\cr
#'   Path to a directory in which Stan should save its output files.
#' @param subset_algos TBD.
#' @param subset_sims TBD.
#' 
#' @returns (`data.table()`)\cr
#'   A data.table containing the simulation results.
#'   
#' @export 
run_sim = function(row_params,
                   num_sims = 1000,
                   constants,
                   num_cores = 4,
                   stan_output_dir = NULL,
                   algo_ids = NULL,
                   rep_ids = NULL) {
  # Sanity checks
  chk$assert_data_table(row_params, nrows = 1)
  chk$assert_subset(
    colnames(row_params),
    c("scenario.id", "n_trt", "n_ctrl", "n_hist", "b_trt", "shift_x2", "shift_x3", "shift_u")
  )
  chk$assert_set_equal(
    setdiff(colnames(row_params), c("n_trt", "n_ctrl")),
    c("scenario.id", "n_hist", "b_trt", "shift_x2", "shift_x3", "shift_u")
  )
  chk$assert_count(num_sims, positive = TRUE)
  chk$assert_list(constants, names = "unique")
  chk$assert_count(num_cores, positive = TRUE)

  
  # If specified/requested, create sub-directory for Stan output
  if (!is.null(stan_output_dir)) {
    chk$assert_directory_exists(stan_output_dir)
    # For each simulation/scenario that we run, we create a sub-directory 
    # that gets deleted later.
    stan_output_subdir = fs$path(
      stan_output_dir,
      row_params$scenario.id
    )
    fs$dir_create(stan_output_subdir)
    on.exit(fs$dir_delete(stan_output_subdir), add = TRUE)
  }
  
  
  # Either run simulation for selected subset of algorithms or all of them
  if (!is.null(algo_ids)) {
    chk$assert_subset(algo_ids, choices = 1:5)
  } else {
    algo_ids = 1:5
  }
  
  
  # Initialize simulator for generating the data sets
  simulator = Simulator$new(li_params = row_params[, -c("scenario.id")])
  
  # Generate the data sets
  li_data = with_seed(row_params$scenario.id, {
    simulator$simulate(num_sims)
  }, .rng_kind = "Mersenne-Twister")
  
  # Subset the data sets if requested (i.e. only run certain repetitions)
  if (!is.null(rep_ids)) {
    chk$assert_subset(rep_ids, choices = seq_len(num_sims))
    li_data = li_data[rep_ids]
  }
  
  
  # Initialize Stan model for PS-commensurate
  if (4 %in% algo_ids) {
    stan_ps_comm = init_commensurate(
      data_init = li_data[[1]],
      outcome = "y", treatment = "trt", source = "source",
      covariates = c("x1", "x2", "x3"),
      family = gaussian()
    )
  }
  
  
  # Apply different algorithms to the generated data sets
  # Note: Results for Stan models will not be *exactly* reproducible when using 
  #       different number of cores and/or running the code on a different computer.
  res_all_reps = with_seed(row_params$scenario.id, {

    mclapply(seq_along(li_data), function(i) {
      res_one_rep = wrapper_all_algos(
        .data = li_data[[i]],
        constants = constants,
        stan_model = stan_ps_comm,
        stan_output_dir = stan_output_subdir,
        algo_ids = algo_ids
      )
      res_one_rep[, rep.id := i]
      return(res_one_rep)
    }, mc.cores = min(c(num_cores, availableCores())))

  }, .rng_kind = "L'Ecuyer-CMRG")
  
  
  # Combine and organize the results
  res_all_reps = rbindlist(res_all_reps)
  res_all_reps[, scenario.id := row_params$scenario.id]
  setcolorder(
    res_all_reps,
    c("scenario.id", "rep.id", "algo.id", "est", "pval", "ci_lower", "ci_upper", "rhat", "borrow_stat")
  )
  
  
  # Done
  return(res_all_reps[])
}
