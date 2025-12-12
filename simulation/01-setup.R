options(box.path = "R")

box::use(
  data.table[...],
  fs,
  bt = batchtools
)

box::use(
  simfuns/run_sim[run_sim],
  simfuns/gen_data[make_scenarios]
)


# Create simulation registry
reg = bt$makeRegistry(
  file.dir = fs$path_home("scratch", "simregistry"),
  work.dir = fs$path_wd(),
  seed = 1
)

# Specify Slurm scheduler backend
reg$cluster.functions = bt$makeClusterFunctionsSlurm(
  template = fs$path("simulation", "slurm", ext = "tmpl"),
  array.jobs = TRUE
)

# Create directory for Stan output
stan_output_dir = fs$path_home("scratch", "stan_output")
fs$dir_create(stan_output_dir)

# Specify all relevant parameters and function arguments
scenarios = make_scenarios()

constants = list(
  ci_level = 0.95,
  num_warmup = 2500,
  num_posterior = 2500,
  num_chains = 4
)

fixed_args = list(
  num_sims = 1000,
  constants = constants,
  num_cores = 1,
  stan_output_dir = stan_output_dir
)

# Specify the jobs to be submitted later
bt$batchMap(
  run_sim,
  row_params = split(scenarios, by = "scenario.id"),
  more.args = fixed_args,
  reg = reg
)

# Clean the environment
rm(list = ls())
box::purge_cache()
