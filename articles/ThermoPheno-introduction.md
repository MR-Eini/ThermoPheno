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
