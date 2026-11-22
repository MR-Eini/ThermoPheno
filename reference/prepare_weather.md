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

## Details

The bundled files ending in `_dummy_data.csv` contain synthetic weather
for demonstrations. See
`system.file("extdata", "DWD_DATA_LICENSE.md", package = "ThermoPheno")`
for the source and license of the separate DWD example.
