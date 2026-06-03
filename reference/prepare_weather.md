# Prepare daily weather data for ThermoPheno

Standardises column names and adds `year`, `doy`, and `tmean` columns.

## Usage

``` r
prepare_weather(df)
```

## Arguments

- df:

  Data frame with at least `date`, `tmin`, and `tmax` columns.

## Value

A data frame sorted by date.
