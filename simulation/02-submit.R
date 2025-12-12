options(box.path = "R")

box::use(
  data.table[...],
  fs,
  bt = batchtools
)


# Load registry
reg = bt$loadRegistry(fs$path_home("scratch", "simregistry"), writeable = TRUE)

# Specify Slurm scheduler backend
reg$cluster.functions = bt$makeClusterFunctionsSlurm(
  template = fs$path("simulation", "slurm", ext = "tmpl"),
  array.jobs = TRUE
)

# Specify jobs to submit and chunk them for using job arrays
# 1116 scenarios -> (roughly) 35 chunks/arrays of size 32
to_submit = bt$getJobTable(reg = reg)[, .(job.id)]
to_submit[, chunk := bt$chunk(job.id, chunk.size = 32, shuffle = FALSE)]

# Specify Slurm resources
slurm_resources = list(
  qos = "1d",
  # Time in seconds
  walltime = 60 * 60 * 11,
  # Memory in Mb
  memory = 3000,
  ncpus = 1,
  chunks.as.arrayjobs = TRUE
)

# Submit the jobs
bt$submitJobs(
  ids = to_submit,
  resources = slurm_resources,
  reg = reg
)

# Clean the environment
rm(list = ls())
box::purge_cache()
