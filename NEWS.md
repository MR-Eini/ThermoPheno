# ThermoPheno 0.1.1

- Documented the DWD validation workflow and its held-out evaluation.
- Added an attributed daily weather example from a DWD station.
- Extended tests for thermal-time boundaries, planting windows, input handling
  and simulation outcomes.

# ThermoPheno 0.1.0

- Initial formal R-package structure.
- Added thermal-time functions: simple, capped, and triangular methods.
- Added summer- and winter-crop simulation logic.
- Added Shiny application launcher.
- Added observed-calendar thermal-time calibration through `estimate_required_tt_from_observed()`.
- Added testthat unit tests and example weather/scenario data.
