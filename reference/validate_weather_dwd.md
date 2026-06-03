# Validate weather input against ThermoPheno and DWD-like conventions

Performs structural checks for weather tables used by ThermoPheno and
adds practical quality checks inspired by typical Deutscher Wetterdienst
daily temperature datasets.

## Usage

``` r
validate_weather_dwd(df, strict = FALSE)
```

## Arguments

- df:

  A data frame containing at least `date`, `tmin`, and `tmax`.

- strict:

  If `TRUE`, stop on validation errors. If `FALSE`, return issues.

## Value

A list with `ok`, `errors`, and `warnings`.
