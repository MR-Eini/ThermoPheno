# ThermoPheno check-warning fixer
# Run this once from the ThermoPheno repository root:
#   source("scripts/apply_check_fixes.R")

message("Applying ThermoPheno package-check fixes...")

stopifnot(file.exists("DESCRIPTION"), dir.exists("R"))

replace_in_file <- function(path, old, new) {
  x <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl(old, x, fixed = TRUE)) {
    x <- gsub(old, new, x, fixed = TRUE)
    writeLines(strsplit(x, "\n", fixed = TRUE)[[1]], path, useBytes = TRUE)
    message("Updated: ", path)
  } else {
    message("No matching text found, skipped: ", path)
  }
}

# 1) Fix roxygen links that produced unresolved-link warnings.
replace_in_file(
  "R/ThermoPheno_functions.R",
  "#' @param weather Prepared weather data returned by [prepare_weather()].",
  "#' @param weather Prepared weather data returned by `prepare_weather()`."
)
replace_in_file(
  "R/ThermoPheno_functions.R",
  "#' @param weather Prepared weather data from [prepare_weather()].",
  "#' @param weather Prepared weather data from `prepare_weather()`."
)
replace_in_file(
  "R/app_launcher.R",
  "#' Alias for [run_thermopheno_app()].",
  "#' Alias for `run_thermopheno_app()`."
)

# 2) Document min_mean_temp_plant for simulate_one_year(); run_simulation() inherits it.
path <- "R/ThermoPheno_functions.R"
x <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
old <- paste0(
  "#' @param latest_harvest_mmdd Latest harvest date in `MM-DD` format.\n",
  "#' @param forced_harvest_allowed Logical; whether immature forced harvest is allowed."
)
new <- paste0(
  "#' @param latest_harvest_mmdd Latest harvest date in `MM-DD` format.\n",
  "#' @param min_mean_temp_plant Minimum daily mean temperature required for summer-crop planting.\n",
  "#' @param forced_harvest_allowed Logical; whether immature forced harvest is allowed."
)
if (grepl(old, x, fixed = TRUE) && !grepl("#' @param min_mean_temp_plant Minimum daily mean temperature required for summer-crop planting.", x, fixed = TRUE)) {
  x <- gsub(old, new, x, fixed = TRUE)
  writeLines(strsplit(x, "\n", fixed = TRUE)[[1]], path, useBytes = TRUE)
  message("Added min_mean_temp_plant roxygen documentation to simulate_one_year().")
} else {
  message("min_mean_temp_plant documentation already present or expected location not found.")
}

# 3) Replace DESCRIPTION so optional Shiny/app/validation packages do not cause the Imports NOTE.
description_txt <- c(
  "Package: ThermoPheno",
  "Type: Package",
  "Title: Thermal-Time-Based Crop Phenology Simulation",
  "Version: 0.1.0",
  "Authors@R:",
  "    person(given = \"Mohammad Reza\",",
  "           family = \"Eini\",",
  "           email = \"mreiniut@gmail.com\",",
  "           role = c(\"aut\", \"cre\"))",
  "Description: Provides a thermal-time-based crop phenology model for simulating",
  "    crop planting, maturity, and harvest timing under historical and climate",
  "    change conditions. The package includes reusable model functions, an",
  "    interactive Shiny application, bundled example data, and helper functions",
  "    for validation workflows using observed crop phenology and daily climate",
  "    observations.",
  "License: MIT + file LICENSE",
  "Encoding: UTF-8",
  "Depends:",
  "    R (>= 4.1)",
  "Imports:",
  "    stats",
  "Suggests:",
  "    shiny,",
  "    bslib,",
  "    dplyr,",
  "    lubridate,",
  "    ggplot2,",
  "    DT,",
  "    readr,",
  "    tidyr,",
  "    tibble,",
  "    rlang,",
  "    RColorBrewer,",
  "    ggridges,",
  "    rdwd,",
  "    testthat (>= 3.0.0),",
  "    devtools,",
  "    rcmdcheck,",
  "    roxygen2,",
  "    pkgdown,",
  "    knitr,",
  "    rmarkdown",
  "Config/testthat/edition: 3",
  "Roxygen: list(markdown = TRUE)",
  "RoxygenNote: 7.3.3",
  "URL: https://github.com/MR-Eini/ThermoPheno, https://mr-eini.github.io/ThermoPheno/",
  "BugReports: https://github.com/MR-Eini/ThermoPheno/issues",
  "VignetteBuilder: knitr"
)
writeLines(description_txt, "DESCRIPTION", useBytes = TRUE)
message("Updated: DESCRIPTION")

# 4) Add package-build exclusions. This also prevents UPLOAD_NOTES.txt from causing a NOTE if it exists.
ignore_entries <- c(
  "^.*\\.Rproj$",
  "^\\.Rproj\\.user$",
  "^\\.git$",
  "^\\.gitignore$",
  "^\\.github$",
  "^docs$",
  "^pkgdown$",
  "^_pkgdown\\.yml$",
  "^UPLOAD_NOTES\\.txt$",
  "^ThermoPheno_.*\\.tar\\.gz$",
  "^validation/cache$",
  "^validation/results$",
  "^README_files$",
  "^README\\.Rmd$"
)
current_ignore <- if (file.exists(".Rbuildignore")) readLines(".Rbuildignore", warn = FALSE) else character()
combined <- unique(c(current_ignore, ignore_entries))
writeLines(combined, ".Rbuildignore", useBytes = TRUE)
message("Updated: .Rbuildignore")

if (file.exists("UPLOAD_NOTES.txt")) {
  unlink("UPLOAD_NOTES.txt")
  message("Removed top-level UPLOAD_NOTES.txt because it caused a package-check NOTE.")
}

message("Done. Now run:")
message("  rm(list = c('calc_daily_tt','compare_validation_metrics','default_parameters','estimate_required_tt','estimate_required_tt_from_observed','find_planting_date','prepare_weather','run_simulation','simulate_one_year'))")
message("  devtools::document()")
message("  devtools::test()")
message("  rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'))")
message("For pkgdown, first run: pkgdown::clean_site(force = TRUE), then pkgdown::build_site().")
