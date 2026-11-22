# Synthetic calibrate-then-simulate round trip, runnable against any installed
# ThermoPheno version. With 0.1.1 the only calibration available is the raw sum;
# with >= 0.2.0 both accumulation settings are run.
# Rscript validation/diagnostics/roundtrip/roundtrip_synthetic.R
suppressPackageStartupMessages(library(ThermoPheno))
prepare_weather <- ThermoPheno::prepare_weather
simulate_one_year <- ThermoPheno::simulate_one_year
calc_daily_tt <- ThermoPheno::calc_daily_tt
source("tests/testthat/helper-synthetic-weather.R")
weather <- synthetic_seasonal_weather()
cases <- list(
  winter = list(rules = winter_test_rules, plant = as.Date("2000-10-10"),
                harvest = as.Date("2001-07-25"), year = 2001L, deadline = "08-15"),
  summer = list(rules = summer_test_rules, plant = as.Date("2001-05-01"),
                harvest = as.Date("2001-09-20"), year = 2001L, deadline = "11-15")
)
has_acc <- "accumulation" %in% names(formals(estimate_required_tt_from_observed))
modes <- if (has_acc) c("raw", "simulation") else "raw (0.1.1 behaviour)"
out <- list()
for (nm in names(cases)) for (m in modes) {
  cs <- cases[[nm]]
  obs <- data.frame(crop_year = cs$year, observed_planting_date = cs$plant,
                    observed_harvest_date = cs$harvest)
  cal_args <- list(weather = weather, observed_calendar = obs, calibration_years = cs$year,
                   t_base = cs$rules$t_base, t_opt = cs$rules$t_opt,
                   t_max_cut = cs$rules$t_max_cut, tt_mode = cs$rules$tt_mode)
  if (has_acc) cal_args <- c(cal_args, cs$rules[setdiff(names(cs$rules), names(cal_args))],
                             list(accumulation = m))
  req <- do.call(estimate_required_tt_from_observed, cal_args)$required_tt
  sim <- simulate_from_planting(weather, cs$rules, cs$plant, cs$year, req, cs$deadline)
  out[[paste(nm, m)]] <- data.frame(
    version = as.character(packageVersion("ThermoPheno")), crop = nm, calibration = m,
    required_tt = sprintf("%.10f", req), observed_harvest = as.character(cs$harvest),
    simulated_status = sim$status, simulated_maturity = as.character(sim$maturity_date),
    maturity_error_days = as.numeric(sim$maturity_date - cs$harvest),
    accumulated_tt_at_deadline_or_maturity = sim$accumulated_tt,
    round_trip_within_1_day = isTRUE(abs(as.numeric(sim$maturity_date - cs$harvest)) <= 1))
}
res <- do.call(rbind, out); rownames(res) <- NULL
options(width = 200); print(res)
