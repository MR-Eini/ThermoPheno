# Estimate required thermal time from observed planting and harvest dates

For each observed season, thermal time is accumulated from the observed
planting date to the observed harvest date (both inclusive). The
required thermal time is the median (or mean) of these seasonal sums.

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
  min_weather_coverage = 0.95,
  crop_type = c("summer", "winter"),
  winter_dormancy_temp = 0,
  vernalization_required = FALSE,
  vernalization_temp_min = 0,
  vernalization_temp_max = 10,
  vernalization_days_required = 30,
  spring_regrowth_temp = 5,
  accumulation = c("simulation", "raw")
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

- crop_type:

  `summer` or `winter`. Selects the accumulation rules used with
  `accumulation = "simulation"`; pass the same value as to the
  simulation.

- winter_dormancy_temp, vernalization_required, vernalization_temp_min,
  vernalization_temp_max, vernalization_days_required,
  spring_regrowth_temp:

  Winter-crop rules, as in
  [`simulate_one_year()`](https://mr-eini.github.io/ThermoPheno/reference/simulate_one_year.md).
  Pass the values used for the simulation. Ignored for summer crops and
  with `accumulation = "raw"`.

- accumulation:

  `simulation` (default) accumulates thermal time with the simulation
  rules; `raw` reproduces the 0.1.1 raw sum. See the Accumulation
  section.

## Value

A list containing yearly observed TT values and the estimated required
TT, plus the `accumulation` and `crop_type` used.

## Accumulation

With `accumulation = "simulation"` (the default), each season is
accumulated with the same code and rules that
[`simulate_one_year()`](https://mr-eini.github.io/ThermoPheno/reference/simulate_one_year.md)
uses after planting: the thermal-time method and, for
`crop_type = "winter"`, the vernalization and dormancy rules. A
requirement calibrated this way is reached by the simulation on the
observed harvest date when the simulation is started on the observed
planting date of the same season.

With `accumulation = "raw"`, the daily thermal time is summed over the
observed interval without the winter-crop rules. This reproduces
ThermoPheno 0.1.1. For winter crops with vernalization, the raw sum is
larger than what the simulation accumulates over the same interval, so a
raw requirement is systematically too high for the simulation to reach
on time.

For summer crops (the default `crop_type`), both settings sum the same
daily values, so calls written for ThermoPheno 0.1.1 return the same
values up to floating-point rounding (the simulation adds day by day in
double precision; the raw sum uses
[`sum()`](https://rdrr.io/r/base/sum.html)). That rounding difference
could, in 0.1.1, place a simulated maturity one day after the observed
harvest date of a calibration season.
