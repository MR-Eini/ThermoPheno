# Changelog

## ThermoPheno 0.2.0

### Calibration consistent with the simulation

- [`estimate_required_tt_from_observed()`](https://mr-eini.github.io/ThermoPheno/reference/estimate_required_tt_from_observed.md)
  now accumulates each observed season with the same code and rules as
  [`simulate_one_year()`](https://mr-eini.github.io/ThermoPheno/reference/simulate_one_year.md),
  including the vernalization and dormancy rules for winter crops. New
  arguments: `crop_type`, `winter_dormancy_temp`,
  `vernalization_required`, `vernalization_temp_min`,
  `vernalization_temp_max`, `vernalization_days_required`,
  `spring_regrowth_temp` and `accumulation = c("simulation", "raw")`.
  The returned list also reports `accumulation` and `crop_type`.
- In 0.1.1 the observed-calendar calibration summed raw thermal time.
  For winter crops with vernalization, that sum is larger than what the
  simulation accumulates over the same season, so the calibrated
  requirement was too high and simulated harvests were late or forced at
  the deadline. Calibrating a winter crop with `crop_type = "winter"`
  and the simulation’s rules corrects this. `accumulation = "raw"`
  reproduces the 0.1.1 values.
- For summer crops, and for calls that do not set `crop_type`, results
  equal 0.1.1 up to floating-point rounding (about 1e-12 °C·day),
  because the calibration now adds daily values in the same order as the
  simulation.
- [`estimate_required_tt()`](https://mr-eini.github.io/ThermoPheno/reference/estimate_required_tt.md),
  [`simulate_one_year()`](https://mr-eini.github.io/ThermoPheno/reference/simulate_one_year.md)
  and
  [`run_simulation()`](https://mr-eini.github.io/ThermoPheno/reference/run_simulation.md)
  use one internal accumulator. Simulation results are unchanged;
  summer-crop
  [`estimate_required_tt()`](https://mr-eini.github.io/ThermoPheno/reference/estimate_required_tt.md)
  values can differ from 0.1.1 by floating-point rounding.

### Validation

- Reran the DWD validation with the same protocol and inputs. Winter
  wheat is now calibrated with its winter rules. Forced winter wheat
  harvests fell from 174 to 22 and failed seasons from 8 to 0. Winter
  wheat harvest MAE changed from 14.52 to 9.45 days. The 0.1.1 outputs
  are archived and compared in the validation article.
- Added a calibration-period diagnostic
  (`validation/diagnose_tt_consistency.R`) and run comparison scripts.
- Added round-trip, equality and raw-mode tests for calibration and
  simulation.

## ThermoPheno 0.1.1

- Documented the DWD validation workflow and its held-out evaluation.
- Added an attributed daily weather example from a DWD station.
- Extended tests for thermal-time boundaries, planting windows, input
  handling and simulation outcomes.

## ThermoPheno 0.1.0

- Initial formal R-package structure.
- Added thermal-time functions: simple, capped, and triangular methods.
- Added summer- and winter-crop simulation logic.
- Added Shiny application launcher.
- Added observed-calendar thermal-time calibration through
  [`estimate_required_tt_from_observed()`](https://mr-eini.github.io/ThermoPheno/reference/estimate_required_tt_from_observed.md).
- Added testthat unit tests and example weather/scenario data.
