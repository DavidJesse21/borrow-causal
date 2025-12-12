box::use(
  fs,
  knitr[write_bib]
)

pkgs = c(
  "base",
  "cmdstanr",
  "psborrow2",
  "dbarts",
  "renv",
  "batchtools"
)

file = fs$path("simulation", "packages", ext = "bib")

write_bib(pkgs, file)
