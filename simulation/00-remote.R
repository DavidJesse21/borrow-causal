if (!nzchar(system.file(package = "renv"))) {
  message("Installing `renv` package...")
  install.packages("renv")
}

message("Restoring `renv` project...")
to_install = readRDS("simulation/packages.rds")
renv::restore(packages = to_install)

message("Updating lockfile...")
renv::snapshot()

message("Disabling automatic snapshots and ignoring testthat package...")
options(renv.config.auto.snapshot = FALSE)
rprofile = 'options(
  renv.config.auto.snapshot = FALSE,
  renv.settings.ignored.packages = "testthat"
)'
rprofile = c(rprofile, readLines(".Rprofile"))
writeLines(rprofile, con = ".Rprofile")

cmdstan_installed = tryCatch(
  nzchar(cmdstanr::cmdstan_path()),
  err = \(e) FALSE
)
if (cmdstan_installed) {
  message("cmdstan already installed!")
} else {
  message("Installing cmdstan...")
  cmdstanr::install_cmdstan(version = "2.36.0")
}

message("All done!")
