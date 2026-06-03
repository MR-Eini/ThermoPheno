# Compare observed and simulated phenological dates

Compare observed and simulated phenological dates

## Usage

``` r
compare_validation_metrics(
  df,
  observed_col = "observed_date",
  simulated_col = "simulated_date"
)
```

## Arguments

- df:

  Data frame containing observed and simulated date columns.

- observed_col:

  Name of observed date column.

- simulated_col:

  Name of simulated date column.

## Value

One-row data frame with validation metrics.
