# ThermoPheno

**ThermoPheno** is an R package and Shiny application for thermal-time-based crop phenology simulation. It estimates crop planting, maturity, and harvest timing from daily minimum and maximum temperature data and supports historical analysis, climate-scenario comparison, observed-phenology validation, and benchmarking against external crop calendar or agro-hydrological model outputs.

Package website: <https://mr-eini.github.io/ThermoPheno/>  
Repository: <https://github.com/MR-Eini/ThermoPheno>

## Main features

- **Thermal-time simulation** using simple, capped, or triangular temperature-response formulations.
- **Summer and winter crop workflows**, including simplified winter dormancy and vernalization logic.
- **Baseline thermal-time calibration** from a reference planting date and days-to-maturity period.
- **Dynamic planting-date detection** based on temperature conditions within a user-defined planting window.
- **Climate-scenario comparison** for multi-model, multi-scenario, multi-period, and multi-station files.
- **Observed-phenology validation support** for comparing simulated crop calendar dates with independent phenological observations.
- **Interactive Shiny interface** for uploading weather data, configuring crop parameters, running simulations, and visualising results.

## Installation

### Install from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("MR-Eini/ThermoPheno", dependencies = TRUE)
```

### Install from a local source folder

```r
remotes::install_local("ThermoPheno", dependencies = TRUE)
```

### Install from a source archive

```r
install.packages("ThermoPheno_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Launch the Shiny app

```r
library(ThermoPheno)
run_thermopheno_app()
```

The package-level alias also works:

```r
ThermoPheno()
```

## Example data

Bundled example files are available in `inst/extdata`:

```r
system.file("extdata", "Germany_historical_1981_2010_dummy_data.csv", package = "ThermoPheno")
system.file("extdata", "Germany_10_scenarios_2071_2100_dummy_data.csv", package = "ThermoPheno")
```

## Input data format

ThermoPheno expects daily weather data in CSV format with a header row. Column names are treated case-insensitively.

### Required columns

| Column | Type | Description |
|---|---|---|
| `date` | Date | Daily date in `YYYY-MM-DD` format |
| `tmin` | Numeric | Daily minimum air temperature in °C |
| `tmax` | Numeric | Daily maximum air temperature in °C |

### Optional grouping columns for climate-scenario files

| Column | Description |
|---|---|
| `scenario` | Emissions or climate scenario label, for example `SSP2-4.5` |
| `model` | Climate model or ensemble member name |
| `period` | Time-slice label, for example `2041-2070` |
| `station` | Site or station identifier for multi-station files |

When grouping columns are present, ThermoPheno runs simulations separately for each unique grouping combination and pools the results for comparison and visualisation.

## Minimal model example

```r
library(ThermoPheno)

weather_file <- system.file(
  "extdata", "Germany_historical_1981_2010_dummy_data.csv",
  package = "ThermoPheno"
)

weather <- prepare_weather(read.csv(weather_file))

req <- estimate_required_tt(
  weather = weather,
  baseline_years = 1981:2010,
  planting_mmdd = "04-15",
  days_to_maturity = 140,
  t_base = 8,
  t_opt = 25,
  t_max_cut = 35,
  tt_mode = "triangular",
  crop_type = "summer"
)

sim <- run_simulation(
  weather = weather,
  crop_name = "Maize",
  required_tt = req$required_tt,
  earliest_planting_mmdd = "03-15",
  latest_planting_mmdd = "05-31",
  latest_harvest_mmdd = "10-01",
  t_base = 8,
  t_opt = 25,
  t_max_cut = 35,
  tt_mode = "triangular",
  crop_type = "summer",
  min_mean_temp_plant = 8
)

head(sim)
```

## Main exported functions

| Function | Role |
|---|---|
| `prepare_weather()` | Standardises weather data and adds `year`, `doy`, and `tmean`. |
| `calc_daily_tt()` | Calculates daily thermal time using simple, capped, or triangular methods. |
| `estimate_required_tt()` | Estimates crop-specific required thermal time from a reference baseline period. |
| `estimate_required_tt_from_observed()` | Estimates required thermal time using observed phenological dates. |
| `find_planting_date()` | Finds the first temperature-valid planting date within an allowed window. |
| `simulate_one_year()` | Simulates one crop season for a selected crop and year. |
| `run_simulation()` | Runs multi-year crop phenology simulations. |
| `default_parameters()` | Returns default parameter sets for supported crop types. |
| `compare_validation_metrics()` | Summarises validation performance between simulated and observed dates. |
| `validate_weather_dwd()` | Supports validation workflows using DWD weather and phenology data. |
| `run_thermopheno_app()` | Launches the interactive Shiny application. |
| `ThermoPheno()` | Short alias for launching the app. |

## Output variables

Each simulated season can include the following outputs:

| Variable | Description |
|---|---|
| `crop_name` | Crop name used in the simulation. |
| `crop_type` | Summer or winter crop type. |
| `year` | Simulation year or harvest year. |
| `planting_date` | Date on which planting conditions were first met. |
| `maturity_date` | Date on which accumulated thermal time reached the required value. |
| `harvest_date` | Simulated harvest date or forced harvest date. |
| `season_length_days` | Number of calendar days from planting to harvest. |
| `accumulated_tt` | Total accumulated thermal time in °C-days. |
| `required_tt` | Required thermal time used for maturity simulation. |
| `maturity_fraction` | Ratio between accumulated and required thermal time. |
| `status` | Simulation status, such as `mature`, `forced_harvest_immature`, `failed_to_mature`, `insufficient_vernalization`, or `not_planted`. |
| `vernalization_days` | Number of days satisfying the vernalization criterion for winter crops. |
| `vernalization_satisfied` | Whether the required vernalization condition was satisfied. |

## Thermal-time methods

### Simple

```text
TT = max(Tmean - Tbase, 0)
```

### Capped

```text
TT = max(min(Tmean, Topt) - Tbase, 0)
```

### Triangular

```text
TT = 0                                                   if Tmean <= Tbase
TT = Tmean - Tbase                                      if Tbase < Tmean <= Topt
TT = (Topt - Tbase) * (Tmax_cut - Tmean) /              if Topt < Tmean < Tmax_cut
     (Tmax_cut - Topt)
TT = 0                                                   if Tmean >= Tmax_cut
```

## Main assumptions

- Crop development is driven by air temperature only.
- Daily mean temperature is calculated as `(tmin + tmax) / 2`.
- The winter-crop workflow uses simplified vernalization and dormancy rules.
- Photoperiod, radiation, water stress, cultivar differences, management adaptation, and frost-kill processes are not explicitly modelled.
- Each growing season is simulated independently.

## Local checks before upload

```r
devtools::document()
devtools::test()
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))
pkgdown::build_site()
```

## Citation

If you use ThermoPheno in a report, repository, or paper, cite it as:

> Eini, M. R. (2026). **ThermoPheno: Thermal-Time-Based Crop Phenology Model**. R package version 0.1.0. GitHub repository: <https://github.com/MR-Eini/ThermoPheno>

## License

MIT. See the `LICENSE` file for details.
