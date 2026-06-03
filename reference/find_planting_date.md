# Find the first planting date within an allowed planting window

Find the first planting date within an allowed planting window

## Usage

``` r
find_planting_date(
  weather_window,
  earliest_planting_date,
  latest_planting_date,
  crop_type = c("summer", "winter"),
  min_mean_temp_plant = -999,
  winter_plant_temp_min = 5,
  winter_plant_temp_max = 15
)
```

## Arguments

- weather_window:

  Prepared weather data covering the planting window.

- earliest_planting_date:

  Earliest allowed planting date.

- latest_planting_date:

  Latest allowed planting date.

- crop_type:

  `summer` or `winter`.

- min_mean_temp_plant:

  Minimum daily mean temperature for summer-crop planting.

- winter_plant_temp_min:

  Minimum daily mean temperature for winter-crop planting.

- winter_plant_temp_max:

  Maximum daily mean temperature for winter-crop planting.

## Value

A list with `found` and `planting_date`.
