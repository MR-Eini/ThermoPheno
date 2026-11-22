test_that("summer crop simulation reaches maturity under simple constant weather", {
  dates <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "day")
  weather <- prepare_weather(data.frame(date = dates, tmin = 10, tmax = 20))
  sim <- simulate_one_year(
    weather = weather,
    sim_year = 2001,
    required_tt = 100,
    earliest_planting_mmdd = "03-01",
    latest_planting_mmdd = "03-31",
    latest_harvest_mmdd = "10-01",
    t_base = 10,
    tt_mode = "simple",
    crop_type = "summer",
    min_mean_temp_plant = 0,
    crop_name = "Test crop"
  )
  expect_equal(sim$status, "mature")
  expect_equal(as.character(sim$planting_date), "2001-03-01")
  expect_equal(as.character(sim$harvest_date), "2001-03-20")
})

test_that("run_simulation returns one row per year for summer crops", {
  dates <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "day")
  weather <- prepare_weather(data.frame(date = dates, tmin = 10, tmax = 20))
  sim <- run_simulation(
    weather = weather,
    crop_name = "Test crop",
    required_tt = 100,
    earliest_planting_mmdd = "03-01",
    latest_planting_mmdd = "03-31",
    latest_harvest_mmdd = "10-01",
    t_base = 10,
    tt_mode = "simple",
    crop_type = "summer"
  )
  expect_equal(nrow(sim), 3)
})

test_that("planting search includes both date boundaries", {
  dates <- as.Date("2001-03-01") + 0:3
  weather <- data.frame(date = dates, tmin = c(10, 5, 10, 10), tmax = c(10, 5, 10, 10))
  planted <- find_planting_date(weather, dates[2], dates[3], min_mean_temp_plant = 10)
  expect_identical(planted$planting_date, dates[3])
  planted <- find_planting_date(weather, dates[1], dates[3], min_mean_temp_plant = 10)
  expect_identical(planted$planting_date, dates[1])
  no_plant <- find_planting_date(weather, dates[2], dates[2], min_mean_temp_plant = 10)
  expect_false(no_plant$found)
  expect_true(is.na(no_plant$planting_date))
})

test_that("summer simulations distinguish forced harvest and failed seasons", {
  weather <- data.frame(date = as.Date("2001-03-01") + 0:3, tmin = 10, tmax = 10)
  force <- simulate_one_year(weather, 2001, 50, "03-01", "03-01", "03-04",
                             t_base = 0, min_fraction_tt_for_forced_harvest = 0.8)
  expect_identical(force$status, "forced_harvest_immature")
  expect_equal(force$maturity_fraction, 0.8)
  expect_identical(force$harvest_date, as.Date("2001-03-04"))
  failed <- simulate_one_year(weather, 2001, 60, "03-01", "03-01", "03-04",
                              t_base = 0, min_fraction_tt_for_forced_harvest = 0.8)
  expect_identical(failed$status, "failed_to_mature")
  expect_true(is.na(failed$harvest_date))
  dormant <- simulate_one_year(weather, 2001, 50, "03-01", "03-01", "03-04",
                               t_base = 0, min_mean_temp_plant = 11)
  expect_identical(dormant$status, "not_planted")
})

test_that("winter simulations start accumulation when vernalization is satisfied", {
  dates <- c(as.Date("2000-10-01") + 0:1, as.Date("2001-03-01"))
  weather <- data.frame(date = dates, tmin = c(5, 5, 20), tmax = c(5, 5, 20))
  run_winter <- function(required_days) {
    run_simulation(weather, "Wheat", 25, "10-01", "10-01", "03-01", t_base = 0,
                   crop_type = "winter", vernalization_required = TRUE,
                   vernalization_days_required = required_days,
                   forced_harvest_allowed = FALSE)
  }
  failed <- run_winter(3)
  expect_identical(failed$status, "insufficient_vernalization")
  expect_equal(failed$accumulated_tt, 0)
  expect_false(failed$vernalization_satisfied)
  mature <- run_winter(2)
  expect_identical(mature$status, "mature")
  expect_identical(mature$planting_date, dates[1])
  expect_identical(mature$maturity_date, dates[3])
  expect_equal(mature$accumulated_tt, 25)
  expect_true(mature$vernalization_satisfied)
})
