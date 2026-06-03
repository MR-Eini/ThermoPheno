# Estimate required thermal time from a reference baseline period

Estimate required thermal time from a reference baseline period

## Usage

``` r
estimate_required_tt(
  weather,
  baseline_years,
  planting_mmdd,
  days_to_maturity,
  t_base,
  t_opt = NA_real_,
  t_max_cut = NA_real_,
  tt_mode = "simple",
  crop_type = c("summer", "winter"),
  winter_dormancy_temp = 0,
  vernalization_required = FALSE,
  vernalization_temp_min = 0,
  vernalization_temp_max = 10,
  vernalization_days_required = 30,
  spring_regrowth_temp = 5
)
```

## Arguments

- weather:

  Prepared weather data returned by
  [`prepare_weather()`](https://mr-eini.github.io/ThermoPheno/reference/prepare_weather.md).

- baseline_years:

  Integer vector of baseline years.

- planting_mmdd:

  Reference planting date in `MM-DD` format.

- days_to_maturity:

  Reference number of days from planting to maturity.

- t_base:

  Base temperature.

- t_opt:

  Optimum temperature.

- t_max_cut:

  Upper cutoff temperature.

- tt_mode:

  Thermal-time method.

- crop_type:

  `summer` or `winter`.

- winter_dormancy_temp:

  Dormancy threshold for winter crops.

- vernalization_required:

  Logical; whether vernalization is required.

- vernalization_temp_min:

  Minimum temperature for vernalization days.

- vernalization_temp_max:

  Maximum temperature for vernalization days.

- vernalization_days_required:

  Required number of vernalization days.

- spring_regrowth_temp:

  Spring regrowth threshold.

## Value

A list with yearly thermal-time values and the mean required thermal
time.
