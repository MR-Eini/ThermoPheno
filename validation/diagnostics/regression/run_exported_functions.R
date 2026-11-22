# Run the exported model functions on fixed inputs with whichever ThermoPheno
# version is installed and save the outputs, so that two versions can be
# compared with compare_versions.R.
#   Rscript validation/diagnostics/regression/run_exported_functions.R <out.rds>
# Inputs: the bundled example data and three cached DWD climate stations
# (validation/cache), calibration years 1991-2010 only.
suppressPackageStartupMessages(library(ThermoPheno))
out_file <- commandArgs(trailingOnly = TRUE)[[1]]

ext <- function(f) system.file("extdata", f, package = "ThermoPheno")
read_dwd_zip <- function(path) {
  td <- tempfile(); dir.create(td); utils::unzip(path, exdir = td)
  f <- list.files(td, pattern = "^produkt", full.names = TRUE)[1]
  d <- utils::read.table(f, sep = ";", header = TRUE, strip.white = TRUE,
                         na.strings = c("-999", "-9999"))
  d <- d[is.finite(d$TNK) & is.finite(d$TXK), ]
  prepare_weather(data.frame(date = as.Date(as.character(d$MESS_DATUM), "%Y%m%d"),
                             tmin = d$TNK, tmax = d$TXK))
}
dwd_zips <- sort(list.files("validation/cache", pattern = "^tageswerte_KL_.*zip$",
                            full.names = TRUE))
start_year <- as.integer(substr(sub("^tageswerte_KL_[0-9]{5}_", "", basename(dwd_zips)), 1, 4))
dwd_zips <- dwd_zips[start_year <= 1989]
dwd_weather <- list()
for (z in dwd_zips) {
  w <- read_dwd_zip(z); w <- w[w$year >= 1990 & w$year <= 2010, ]
  if (nrow(w) >= 0.95 * 21 * 365) dwd_weather[[basename(z)]] <- w
  if (length(dwd_weather) == 3) break
}
hist_bundled <- utils::read.csv(ext("Germany_historical_1981_2010_dummy_data.csv"))
hist_bundled <- hist_bundled[hist_bundled$station == hist_bundled$station[1], c("date", "tmin", "tmax")]
weathers <- c(
  list(bundled_historical = prepare_weather(hist_bundled)),
  dwd_weather
)

res <- list(version = as.character(packageVersion("ThermoPheno")))
for (wn in names(weathers)) {
  w <- weathers[[wn]]
  yrs <- intersect(1991:2010, unique(w$year))
  for (ct in c("summer", "winter")) {
    p <- default_parameters(ct)
    key <- paste(wn, ct, sep = "|")
    req <- estimate_required_tt(
      w, yrs, p$baseline_planting_mmdd, p$days_to_maturity, p$t_base, p$t_opt,
      p$t_max_cut, "triangular", ct, p$winter_dormancy_temp, p$vernalization_required,
      0, 10, p$vernalization_days_required, p$spring_regrowth_temp)
    sim <- run_simulation(
      w, p$crop_name, req$required_tt, p$earliest_planting_mmdd, p$latest_planting_mmdd,
      p$latest_harvest_mmdd, p$t_base, p$t_opt, p$t_max_cut, "triangular", ct,
      p$min_mean_temp_plant, p$forced_harvest_allowed, p$min_fraction_tt_for_forced_harvest,
      p$winter_dormancy_temp, p$vernalization_required, 0, 10,
      p$vernalization_days_required, p$spring_regrowth_temp,
      p$winter_plant_temp_min, p$winter_plant_temp_max)
    # Observed calendar: the simulated dates, so every function sees the same calendar.
    cal <- data.frame(crop_year = sim$year, observed_planting_date = sim$planting_date,
                      observed_harvest_date = sim$harvest_date)
    # 0.1.1-style call (no crop_type / accumulation arguments).
    obs_req <- tryCatch(estimate_required_tt_from_observed(
      w, cal, yrs, t_base = p$t_base, t_opt = p$t_opt, t_max_cut = p$t_max_cut,
      tt_mode = "triangular"), error = function(e) conditionMessage(e))
    res[[key]] <- list(estimate_required_tt = req, run_simulation = sim,
                       estimate_required_tt_from_observed = obs_req)
  }
}
saveRDS(res, out_file)
cat("Saved", out_file, "for ThermoPheno", res$version, "\n")
