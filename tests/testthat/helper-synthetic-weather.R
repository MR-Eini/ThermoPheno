# Smooth synthetic seasonal climate used by the calibration/simulation
# consistency tests. Daily mean temperature follows an annual sine wave with a
# minimum of about -1 degC in mid-January and a maximum of about 19 degC in
# mid-July; the daily range is 10 degC.
synthetic_seasonal_weather <- function(from = "1999-01-01", to = "2003-12-31") {
  dates <- seq(as.Date(from), as.Date(to), by = "day")
  doy <- as.integer(format(dates, "%j"))
  tmean <- 9 + 10 * sin(2 * pi * (doy - 110) / 365)
  prepare_weather(data.frame(date = dates, tmin = tmean - 5, tmax = tmean + 5))
}

# Winter-crop rules used in the tests (the DWD validation winter wheat setup).
winter_test_rules <- list(
  t_base = 0, t_opt = 18, t_max_cut = 30, tt_mode = "triangular",
  crop_type = "winter", winter_dormancy_temp = 0,
  vernalization_required = TRUE, vernalization_days_required = 45,
  spring_regrowth_temp = 5
)

summer_test_rules <- list(
  t_base = 8, t_opt = 25, t_max_cut = 35, tt_mode = "triangular",
  crop_type = "summer"
)

# Simulate one season with planting forced to `planting_date`.
simulate_from_planting <- function(weather, rules, planting_date, harvest_year,
                                   required_tt, latest_harvest_mmdd,
                                   forced_harvest_allowed = FALSE) {
  mmdd <- format(as.Date(planting_date), "%m-%d")
  do.call(simulate_one_year, c(list(
    weather = weather, sim_year = harvest_year, required_tt = required_tt,
    earliest_planting_mmdd = mmdd, latest_planting_mmdd = mmdd,
    latest_harvest_mmdd = latest_harvest_mmdd,
    min_mean_temp_plant = -999, winter_plant_temp_min = -999,
    winter_plant_temp_max = 999, forced_harvest_allowed = forced_harvest_allowed
  ), rules))
}

# Independent re-implementation of the ThermoPheno 0.1.1 simulation loop
# (simulate_one_year() before the shared accumulator), used as a reference.
reference_v011_simulation_tt <- function(sub, rules) {
  r <- utils::modifyList(list(
    t_opt = NA_real_, t_max_cut = NA_real_, tt_mode = "simple",
    crop_type = "summer", winter_dormancy_temp = 0,
    vernalization_required = FALSE, vernalization_temp_min = 0,
    vernalization_temp_max = 10, vernalization_days_required = 30
  ), rules)
  winter <- r$crop_type == "winter"
  vernal_days <- 0
  vern_sat <- !(winter && isTRUE(r$vernalization_required))
  cum_tt <- 0
  cum_vec <- numeric(nrow(sub))
  for (i in seq_len(nrow(sub))) {
    tmean_i <- sub$tmean[i]
    if (winter && isTRUE(r$vernalization_required) && !vern_sat) {
      if (tmean_i >= r$vernalization_temp_min && tmean_i <= r$vernalization_temp_max) {
        vernal_days <- vernal_days + 1
      }
      if (vernal_days >= r$vernalization_days_required) vern_sat <- TRUE
    }
    if (winter && !vern_sat) {
      tt_i <- 0
    } else {
      tt_i <- calc_daily_tt(sub$tmin[i], sub$tmax[i], r$t_base, r$t_opt,
                            r$t_max_cut, mode = r$tt_mode)
      if (winter && tmean_i <= r$winter_dormancy_temp) tt_i <- 0
    }
    cum_tt <- cum_tt + tt_i
    cum_vec[i] <- cum_tt
  }
  cum_vec
}
