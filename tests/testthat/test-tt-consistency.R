# Calibration and simulation must accumulate thermal time identically.

calibrate_one_season <- function(weather, rules, plant, harvest, year, ...) {
  obs <- data.frame(crop_year = year, observed_planting_date = as.Date(plant),
                    observed_harvest_date = as.Date(harvest))
  do.call(estimate_required_tt_from_observed, c(list(
    weather = weather, observed_calendar = obs, calibration_years = year
  ), rules, list(...)))
}

test_that("round trip: winter crop calibrated from a season matures on its observed harvest date", {
  weather <- synthetic_seasonal_weather()
  plant <- as.Date("2000-10-10")
  harvest <- as.Date("2001-07-25")
  req <- calibrate_one_season(weather, winter_test_rules, plant, harvest, 2001)
  expect_identical(req$accumulation, "simulation")
  sim <- simulate_from_planting(weather, winter_test_rules, plant, 2001,
                                req$required_tt, latest_harvest_mmdd = "08-15")
  expect_identical(sim$status, "mature")
  expect_equal(sim$planting_date, plant)
  expect_lte(abs(as.numeric(sim$maturity_date - harvest)), 1)
})

test_that("round trip: summer crop calibrated from a season matures on its observed harvest date", {
  weather <- synthetic_seasonal_weather()
  plant <- as.Date("2001-05-01")
  harvest <- as.Date("2001-09-20")
  req <- calibrate_one_season(weather, summer_test_rules, plant, harvest, 2001)
  sim <- simulate_from_planting(weather, summer_test_rules, plant, 2001,
                                req$required_tt, latest_harvest_mmdd = "11-15")
  expect_identical(sim$status, "mature")
  expect_lte(abs(as.numeric(sim$maturity_date - harvest)), 1)
})

test_that("raw accumulation reproduces the 0.1.1 calibration and fails the winter round trip", {
  weather <- synthetic_seasonal_weather()
  plant <- as.Date("2000-10-10")
  harvest <- as.Date("2001-07-25")
  raw <- calibrate_one_season(weather, winter_test_rules, plant, harvest, 2001,
                              accumulation = "raw")
  # Value produced by ThermoPheno 0.1.1 for the same inputs
  # (validation/diagnostics/roundtrip/roundtrip_v0.1.1.log).
  expect_equal(raw$required_tt, 2048.8223908266, tolerance = 1e-12)
  sub <- weather[weather$date >= plant & weather$date <= harvest, ]
  expect_identical(raw$required_tt, sum(calc_daily_tt(sub$tmin, sub$tmax, 0, 18, 30,
                                                      mode = "triangular")))
  sim <- simulate_from_planting(weather, winter_test_rules, plant, 2001,
                                raw$required_tt, latest_harvest_mmdd = "08-15")
  expect_gt(as.numeric(sim$maturity_date - harvest), 1)

  summer_raw <- calibrate_one_season(weather, summer_test_rules, "2001-05-01",
                                     "2001-09-20", 2001, accumulation = "raw")
  expect_equal(summer_raw$required_tt, 1224.4860296517, tolerance = 1e-12)
})

test_that("summer calibration is unchanged from 0.1.1 up to floating-point rounding", {
  weather <- synthetic_seasonal_weather()
  obs <- data.frame(
    crop_year = 2000:2003,
    observed_planting_date = as.Date(c("2000-04-25", "2001-05-01", "2002-05-06", "2003-04-28")),
    observed_harvest_date = as.Date(c("2000-09-15", "2001-09-20", "2002-10-01", "2003-09-25"))
  )
  args <- list(weather = weather, observed_calendar = obs, calibration_years = 2000:2003,
               t_base = 8, t_opt = 25, t_max_cut = 35, tt_mode = "triangular")
  default_call <- do.call(estimate_required_tt_from_observed, args)
  raw_call <- do.call(estimate_required_tt_from_observed, c(args, accumulation = "raw"))
  expect_equal(default_call$yearly_required_tt, raw_call$yearly_required_tt, tolerance = 1e-12)
  expect_equal(default_call$required_tt, raw_call$required_tt, tolerance = 1e-12)
})

test_that("winter calibration with simulation rules is below the raw sum", {
  weather <- synthetic_seasonal_weather()
  sim_rules <- calibrate_one_season(weather, winter_test_rules, "2000-10-10", "2001-07-25", 2001)
  raw <- calibrate_one_season(weather, winter_test_rules, "2000-10-10", "2001-07-25", 2001,
                              accumulation = "raw")
  expect_lt(sim_rules$required_tt, raw$required_tt)
})

test_that("shared accumulator equals the simulation's accumulation", {
  weather <- synthetic_seasonal_weather()
  for (cs in list(
    list(rules = winter_test_rules, plant = as.Date("2000-10-10"), end = "08-15"),
    list(rules = summer_test_rules, plant = as.Date("2001-05-01"), end = "11-15"),
    list(rules = utils::modifyList(winter_test_rules, list(vernalization_required = FALSE)),
         plant = as.Date("2000-10-10"), end = "08-15")
  )) {
    end_date <- as.Date(paste0("2001-", cs$end))
    sub <- weather[weather$date >= cs$plant & weather$date <= end_date, ]
    acc <- do.call(ThermoPheno:::accumulate_season_tt,
                   c(list(tmin = sub$tmin, tmax = sub$tmax, tmean = sub$tmean), cs$rules))
    sim <- simulate_from_planting(weather, cs$rules, cs$plant, 2001, Inf, cs$end)
    expect_identical(sim$accumulated_tt, acc$total_tt)
    # The refactor keeps the 0.1.1 simulation arithmetic.
    expect_identical(acc$cum_tt, reference_v011_simulation_tt(sub, cs$rules))
    # Calibration over the same interval uses the same accumulation.
    cal <- calibrate_one_season(weather, cs$rules, cs$plant, end_date, 2001)
    expect_identical(cal$required_tt, acc$total_tt)
  }
})

test_that("baseline calibration uses the shared accumulator", {
  weather <- synthetic_seasonal_weather()
  base <- do.call(estimate_required_tt, c(list(
    weather = weather, baseline_years = 2000:2002, planting_mmdd = "10-10",
    days_to_maturity = 280), winter_test_rules))
  for (yr in 2000:2002) {
    start <- as.Date(paste0(yr, "-10-10"))
    sub <- weather[weather$date >= start & weather$date <= start + 279, ]
    expect_identical(base$yearly_required_tt$required_tt[base$yearly_required_tt$year == yr],
                     reference_v011_simulation_tt(sub, winter_test_rules)[280])
  }
})

test_that("accumulation argument is validated", {
  weather <- synthetic_seasonal_weather()
  expect_error(calibrate_one_season(weather, summer_test_rules, "2001-05-01", "2001-09-20",
                                    2001, accumulation = "other"))
})
