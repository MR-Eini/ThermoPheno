# ThermoPheno introduction

The bundled `*_dummy_data.csv` files contain synthetic weather for
demonstrating the workflow. For an observed weather example, use the
2005 daily temperatures from DWD station 04709 (Simmern-Wahlbach).
Source: Deutscher Wetterdienst (DWD), Climate Data Center, [CC BY
4.0](https://creativecommons.org/licenses/by/4.0/). The extract selects
one year and renames the date and temperature fields.

``` r

library(ThermoPheno)
path <- system.file("extdata", "dwd_simmern_wahlbach_2005.csv", package = "ThermoPheno")
head(prepare_weather(read.csv(path, colClasses = c(station = "character"))))
```

``` r

library(ThermoPheno)
weather <- prepare_weather(read.csv(system.file("extdata", "Germany_historical_1981_2010_dummy_data.csv", package = "ThermoPheno")))
params <- default_parameters("summer")
req <- estimate_required_tt(weather, 1981:2010, params$baseline_planting_mmdd, params$days_to_maturity, params$t_base, params$t_opt, params$t_max_cut, "triangular", "summer")
run_simulation(weather, params$crop_name, req$required_tt, params$earliest_planting_mmdd, params$latest_planting_mmdd, params$latest_harvest_mmdd, params$t_base, params$t_opt, params$t_max_cut, "triangular", "summer")
```

## Calibrating from observed dates

[`estimate_required_tt_from_observed()`](https://mr-eini.github.io/ThermoPheno/reference/estimate_required_tt_from_observed.md)
derives the requirement from observed planting and harvest dates. By
default (`accumulation = "simulation"`), each season is accumulated with
the same rules the simulation uses after planting. Pass the same crop
type and winter-crop rules as to the simulation, so that a simulation
started on an observed planting date matures on the observed harvest
date of that season.

``` r

wheat <- default_parameters("winter")
observed <- data.frame(
  crop_year = 1990:1992,
  observed_planting_date = as.Date(c("1989-10-05", "1990-10-08", "1991-10-03")),
  observed_harvest_date = as.Date(c("1990-08-01", "1991-08-05", "1992-07-28"))
)
req_obs <- estimate_required_tt_from_observed(
  weather, observed, calibration_years = 1990:1992,
  t_base = wheat$t_base, t_opt = wheat$t_opt, t_max_cut = wheat$t_max_cut,
  tt_mode = "triangular", crop_type = "winter",
  winter_dormancy_temp = wheat$winter_dormancy_temp,
  vernalization_required = wheat$vernalization_required,
  vernalization_days_required = wheat$vernalization_days_required,
  spring_regrowth_temp = wheat$spring_regrowth_temp
)
req_obs$required_tt
```

`accumulation = "raw"` sums daily thermal time without the winter-crop
rules, as ThermoPheno 0.1.1 did. For winter crops with vernalization, a
raw requirement is larger than the thermal time the simulation can
accumulate over the same season, which makes simulated harvests late or
forced at the deadline. For summer crops, both settings give the same
values up to floating-point rounding.
