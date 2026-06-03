# Calculate daily thermal time

Calculate daily thermal time

## Usage

``` r
calc_daily_tt(
  tmin,
  tmax,
  t_base,
  t_opt = NA_real_,
  t_max_cut = NA_real_,
  mode = c("simple", "capped", "triangular")
)
```

## Arguments

- tmin:

  Daily minimum temperature in degrees Celsius.

- tmax:

  Daily maximum temperature in degrees Celsius.

- t_base:

  Base temperature.

- t_opt:

  Optimum temperature; required for `capped` and `triangular` modes.

- t_max_cut:

  Upper cutoff temperature; required for `triangular` mode.

- mode:

  One of `simple`, `capped`, or `triangular`.

## Value

Numeric vector of daily thermal time values.
