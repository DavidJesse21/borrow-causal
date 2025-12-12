options(box.path = "R")

box::use(
  bt = batchtools,
  fs,
  data.table[...]
)

reg = bt$loadRegistry(fs$path_home("scratch", "simregistry"))

res = rbindlist(bt$reduceResultsList())

res_dir = fs$path("simulation", "results")
if (!fs$dir_exists(res_dir)) fs$dir_create(res_dir)
saveRDS(res, fs$path(res_dir, paste0("simres_", Sys.Date()), ext = "rds"))
