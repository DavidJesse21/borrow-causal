# Identify those packages that are required for the simulation study.

box::use(
  fs,
  data.table[rbindlist]
)

paths = list(
  algos = fs$path("R", "borrow"),
  sim1 = fs$path("R", "simfuns"),
  sim2 = fs$path("simulation")
)

deps = rbindlist(lapply(paths, renv::dependencies))
deps = unique(deps$Package)
deps = deps[!deps == "testthat"]

saveRDS(deps, fs$path("simulation", "packages.rds"))
