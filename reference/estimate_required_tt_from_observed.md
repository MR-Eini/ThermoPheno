# Estimate required thermal time from observed planting and harvest dates

Estimate required thermal time from observed planting and harvest dates

## Usage

``` r
estimate_required_tt_from_observed(
  weather,
  observed_calendar,
  calibration_years,
  planting_col = "observed_planting_date",
  harvest_col = "observed_harvest_date",
  year_col = "crop_year",
  t_base,
  t_opt = NA_real_,
  t_max_cut = NA_real_,
  tt_mode = "simple",
  summary_fun = c("median", "mean"),
  min_weather_coverage = 0.95
)
```

## Arguments

- weather:

  Prepared weather data from
  [`prepare_weather()`](https://mr-eini.github.io/ThermoPheno/reference/prepare_weather.md).

- observed_calendar:

  Data frame with observed planting and harvest dates.

- calibration_years:

  Years used for thermal-time calibration.

- planting_col:

  Name of observed planting-date column.

- harvest_col:

  Name of observed harvest-date column.

- year_col:

  Name of crop-year column.

- t_base:

  Base temperature.

- t_opt:

  Optimum temperature.

- t_max_cut:

  Upper temperature cutoff.

- tt_mode:

  Thermal-time method: `simple`, `capped`, or `triangular`.

- summary_fun:

  Summary method for annual TT values: `median` or `mean`.

- min_weather_coverage:

  Minimum fraction of daily weather records required.

## Value

A list containing yearly observed TT values and the estimated required
TT.
