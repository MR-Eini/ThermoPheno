# Run ThermoPheno simulations across all years in a weather data set

Run ThermoPheno simulations across all years in a weather data set

## Usage

``` r
run_simulation(
  weather,
  crop_name,
  required_tt,
  earliest_planting_mmdd,
  latest_planting_mmdd,
  latest_harvest_mmdd,
  t_base,
  t_opt = NA_real_,
  t_max_cut = NA_real_,
  tt_mode = "simple",
  crop_type = c("summer", "winter"),
  min_mean_temp_plant = 0,
  forced_harvest_allowed = TRUE,
  min_fraction_tt_for_forced_harvest = 0,
  winter_dormancy_temp = 0,
  vernalization_required = FALSE,
  vernalization_temp_min = 0,
  vernalization_temp_max = 10,
  vernalization_days_required = 30,
  spring_regrowth_temp = 5,
  winter_plant_temp_min = 5,
  winter_plant_temp_max = 15
)
```

## Arguments

- weather:

  Prepared weather data returned by
  [`prepare_weather()`](https://mr-eini.github.io/ThermoPheno/reference/prepare_weather.md).

- crop_name:

  Crop label.

- required_tt:

  Required thermal time.

- earliest_planting_mmdd:

  Earliest planting date in `MM-DD` format.

- latest_planting_mmdd:

  Latest planting date in `MM-DD` format.

- latest_harvest_mmdd:

  Latest harvest date in `MM-DD` format.

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

- forced_harvest_allowed:

  Logical; whether immature forced harvest is allowed.

- min_fraction_tt_for_forced_harvest:

  Minimum maturity fraction for forced harvest.

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

- winter_plant_temp_min:

  Minimum winter-crop planting temperature.

- winter_plant_temp_max:

  Maximum winter-crop planting temperature.

## Value

Data frame with one row per simulated harvest year.
