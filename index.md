# ThermoPheno

**Crop Phenology Simulation Model — Historical Analysis & Climate Change
Impact Assessment** ThermoPheno: A thermal-time-based phenology model
for assessing crop growing season shifts under climate change
ThermoPheno models the timing of crop planting, maturity, and harvest
using daily temperature data and the growing degree day (GDD / thermal
time) framework. It runs as a self-contained R Shiny application and
supports both historical climate analysis and multi-scenario climate
change projections.

------------------------------------------------------------------------

## Features

- **Three thermal time modes** — simple (linear), capped (plateau),
  triangular (heat-stress)
- **Summer and winter crop support** — full vernalization and winter
  dormancy logic for cool-season cereals
- **Automatic calibration** — required thermal time is estimated from a
  user-supplied historical baseline, not lookup tables
- **Dynamic planting date** — planting triggers on temperature
  conditions, not fixed calendar dates
- **Climate scenario comparison** — upload multi-model / multi-scenario
  projections; the app splits and simulates each group
- **Rich visualisations** — ridge/density plots, temperature ribbon
  charts, monthly cycle plots, annual boxplots, season-timing heatmaps

------------------------------------------------------------------------

## Repository Structure

    ThermoPheno/
    ├── DESCRIPTION
    ├── NAMESPACE
    ├── R/
    │   ├── ThermoPheno_functions.R
    │   ├── app_launcher.R
    │   └── zzz.R
    ├── inst/
    │   ├── app/
    │   │   └── app.R
    │   └── extdata/
    │       ├── Germany_historical_1981_2010_dummy_data.csv
    │       └── Germany_10_scenarios_2071_2100_dummy_data.csv
    ├── man/
    │   └── ThermoPheno.Rd
    ├── .github/workflows/
    │   ├── R-CMD-check.yaml
    │   └── pkgdown.yaml
    ├── _pkgdown.yml
    └── README.md

`ThermoPheno_functions.R` and `app.R` share the same modelling logic.
The functions file is designed to be sourced independently for scripted
or batch analyses, while `app.R` is fully self-contained for deployment
as a Shiny app.

------------------------------------------------------------------------

## Quick Start

### Prerequisites

R (≥ 4.1) with the following packages:

``` r
install.packages(c(
  "shiny", "dplyr", "lubridate", "ggplot2",
  "DT", "readr", "ggridges", "RColorBrewer",
  "tidyr", "tibble"
))
```

The app will automatically install any missing packages on first launch.

### Running the App

``` r
shiny::runApp("app.R")
```

Or from within RStudio, open `app.R` and click **Run App**.

------------------------------------------------------------------------

## Input Data Format

Both the historical weather file and the optional climate projection
file are **CSV** with a header row. Column names are case-insensitive.

### Required columns

| Column | Type         | Description                       |
|--------|--------------|-----------------------------------|
| `date` | Date (ISO)   | Daily date in `YYYY-MM-DD` format |
| `tmin` | Numeric (°C) | Daily minimum air temperature     |
| `tmax` | Numeric (°C) | Daily maximum air temperature     |

### Optional grouping columns (climate projection files)

| Column     | Description                                 |
|------------|---------------------------------------------|
| `scenario` | Emissions scenario label (e.g., `SSP2-4.5`) |
| `model`    | Climate model name (e.g., `MPI-ESM`)        |
| `period`   | Time slice label (e.g., `2041-2070`)        |
| `station`  | Site identifier for multi-station files     |

When grouping columns are present, simulations are run separately for
each unique combination and results are pooled for comparison.

------------------------------------------------------------------------

## Pre-Configured Crops

### Maize *(Zea mays)* — Summer crop

| Parameter           | Default        |
|---------------------|----------------|
| Base temperature    | 8 °C           |
| Optimum temperature | 25 °C          |
| Upper cutoff        | 35 °C          |
| Days to maturity    | 140            |
| Reference planting  | 15 April       |
| Planting window     | 1 Mar – 30 Jun |
| Planting threshold  | T_mean ≥ 8 °C  |
| Thermal time mode   | Triangular     |
| Vernalization       | Not required   |

### Winter Wheat *(Triticum aestivum)* — Winter crop

| Parameter               | Default        |
|-------------------------|----------------|
| Base temperature        | 0 °C           |
| Optimum temperature     | 18 °C          |
| Upper cutoff            | 30 °C          |
| Days to maturity        | 300            |
| Reference planting      | 1 October      |
| Planting window         | 1 Sep – 30 Nov |
| Planting temp. window   | 5 – 15 °C      |
| Thermal time mode       | Triangular     |
| Vernalization window    | 0 – 10 °C      |
| Vernalization days req. | 30 days        |
| Dormancy threshold      | ≤ 0 °C         |
| Spring regrowth trigger | 5 °C           |

------------------------------------------------------------------------

## Thermal Time Methods

### Simple

    TT = max(T_mean - T_base, 0)

### Capped

    TT = max(min(T_mean, T_opt) - T_base, 0)

### Triangular

    TT = 0                                              if T_mean ≤ T_base
    TT = T_mean - T_base                               if T_base < T_mean ≤ T_opt
    TT = (T_opt - T_base) × (T_max_cut - T_mean)      if T_opt < T_mean < T_max_cut
             / (T_max_cut - T_opt)
    TT = 0                                              if T_mean ≥ T_max_cut

------------------------------------------------------------------------

## Output Variables

Each simulated season produces:

| Variable                  | Description                                                                                              |
|---------------------------|----------------------------------------------------------------------------------------------------------|
| `planting_date`           | Date planting conditions were first met                                                                  |
| `maturity_date`           | Date accumulated TT ≥ TT_req; `NA` if crop did not mature                                                |
| `harvest_date`            | Actual harvest date (= maturity date, forced date, or `NA`)                                              |
| `season_length_days`      | Calendar days from planting to harvest                                                                   |
| `accumulated_tt`          | Total thermal time accumulated (°C-days)                                                                 |
| `required_tt`             | Calibration-estimated TT for maturity (°C-days)                                                          |
| `maturity_fraction`       | `accumulated_tt / required_tt` (1.0 = fully mature)                                                      |
| `status`                  | `mature` / `forced_harvest_immature` / `failed_to_mature` / `insufficient_vernalization` / `not_planted` |
| `vernalization_days`      | (Winter crops) Days that met the vernalization criterion                                                 |
| `vernalization_satisfied` | (Winter crops) Whether V_days ≥ V_req was reached                                                        |

------------------------------------------------------------------------

## Scripted / Batch Use

Source the functions file and call `run_simulation()` directly:

``` r
source("ThermoPheno_functions.R")

weather <- prepare_weather(read.csv("my_weather.csv"))

# Calibrate
cal <- estimate_required_tt(
  weather          = weather,
  baseline_years   = 1990:2010,
  planting_mmdd    = "04-15",
  days_to_maturity = 140,
  t_base           = 8,
  t_opt            = 25,
  t_max_cut        = 35,
  tt_mode          = "triangular",
  crop_type        = "summer"
)

# Simulate all years
results <- run_simulation(
  weather              = weather,
  crop_name            = "Maize",
  required_tt          = cal$required_tt,
  earliest_planting_mmdd = "03-15",
  latest_planting_mmdd   = "05-31",
  latest_harvest_mmdd    = "10-01",
  t_base               = 8,
  t_opt                = 25,
  t_max_cut            = 35,
  tt_mode              = "triangular",
  crop_type            = "summer"
)

print(results)
```

------------------------------------------------------------------------

## Key Assumptions

- Development is driven by air temperature only (no photoperiod,
  radiation, or water stress).
- Daily mean temperature is the arithmetic mean of T_min and T_max.
- Vernalization is a simple counter — de-vernalization by warm spells is
  not modelled.
- Frost-kill logic is not included.
- Each growing season is simulated independently (no carry-over between
  years).

------------------------------------------------------------------------

## Theoretical Background

The model is grounded in the classical thermal time / growing degree day
literature:

- McMaster & Wilhelm (1997). *Growing degree-days: one equation, two
  interpretations.* Agricultural and Forest Meteorology, 87(4), 291–300.
- Wang & Engel (1998). *Simulation of phenological development of wheat
  crops.* Agricultural Systems, 58(1), 1–24.
- Bonhomme (2000). *Bases and limits to using ‘degree.day’ units.*
  European Journal of Agronomy, 13(1), 1–10.

For full theoretical documentation including all equations, crop
parameter rationale, and model assumptions, see
**`ThermoPheno_Theory.docx`**.

------------------------------------------------------------------------

## License

MIT — free to use, modify, and redistribute with attribution.

## Package usage

``` r
# install.packages("remotes")
remotes::install_github("yourusername/ThermoPheno")
library(ThermoPheno)
ThermoPheno()
```

Example files are bundled in `inst/extdata` and can be accessed with:

``` r
system.file("extdata", package = "ThermoPheno")
```
