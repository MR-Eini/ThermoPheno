# ThermoPheno 0.2.0 – update report

Prepared 22 November 2026. Covers the calibration/simulation thermal-time
mismatch, the fix, the rerun of the DWD validation and its comparison with 0.1.1.
Every number below comes from a run whose script and log are in the repository
(paths given). **Status: the work is in the working tree of `main`, not
committed, not pushed, no pull request and no release yet** (see section 7).

## 1. Diagnosis

Script: `validation/diagnose_tt_consistency.R` (calibration harvest years
1991–2010 only; it re-uses the input/parsing code of the validation script
unchanged and the same station-selection rule). Run against the installed
`v0.1.1` tag. Log and tables: `validation/diagnostics/v0.1.1/`.

Fidelity check: the diagnostic reproduces the published 0.1.1 calibration
exactly (60/60 station requirements and 952 seasonal sums identical, maximum
difference 0).

(a) = raw triangular sum over the observed sowing–harvest interval, as
calibration did in 0.1.1. (b) = thermal time accumulated by
`simulate_one_year()` itself over the same interval, with planting forced to the
observed date and the configured winter rules (vernalization 45 days at
0–10 °C, dormancy at ≤ 0 °C).

| Crop | Seasons | (b)/(a) min | p10 | median | mean | p90 | max | Share < 0.99 | Median (a) − (b), °C·day |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Maize | 279 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 | 0 |
| Sugar beet | 445 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 | 0 |
| Winter wheat | 323 | 0.685 | 0.773 | 0.844 | 0.842 | 0.902 | 0.930 | 100% | 360.8 |

Per-station medians for winter wheat range from 0.742 to 0.884 (22 stations,
`tt_consistency_by_station.csv`). One maize and one wheat season were excluded
because the sowing year did not match the crop-year convention.

Reachability before the harvest deadline (simulation's own planting rules,
calibration years): with the 0.1.1 (raw) requirement, winter wheat could not
accumulate the requirement by 15 August in **85.1%** of 343 station-years; with a
requirement calibrated by simulation rules, 23.9%. Maize 10.9% and sugar beet
3.9% under both.

**Conclusion: the mismatch is confirmed for winter wheat only.** For maize and
sugar beet, calibration and simulation already agreed (up to floating-point
rounding; see section 3).

Harvest deadline check (calibration years only, the 25 selected stations):
11.0% of 500 observed winter wheat harvests were after 15 August (median 4 Aug,
p90 16 Aug, p95 20 Aug, latest 4 Sep). I did not judge this clearly
inconsistent with German practice and **left the deadline unchanged** (open
question 1).

## 2. The fix

- New internal function `accumulate_season_tt()` (`R/ThermoPheno_functions.R`):
  the single implementation of daily thermal time + vernalization + dormancy +
  regrowth state, summed sequentially in double precision exactly as the 0.1.1
  simulation loop did.
- `simulate_one_year()` (and therefore `run_simulation()`) and
  `estimate_required_tt()` now call it instead of their own loops.
- `estimate_required_tt_from_observed()` gains `crop_type`,
  `winter_dormancy_temp`, `vernalization_required`, `vernalization_temp_min`,
  `vernalization_temp_max`, `vernalization_days_required`,
  `spring_regrowth_temp` and `accumulation = c("simulation", "raw")`. The
  default is `"simulation"`; `"raw"` keeps the 0.1.1 sum. The returned list also
  has `accumulation` and `crop_type`. Help page, `NEWS.md` and the introduction
  article document it.
- `validation/dwd_validate_thermopheno.R`: calibration and simulation now take
  the thermal-time rules from one list (`tt_rules`), and calibration uses
  `accumulation = "simulation"`. Stations, years, filters, matching,
  minimum-data rules, metrics and `crop_config.csv` are unchanged. The column
  `tt_calibration_method` now reads `observed_calendar/simulation`.

**Version 0.2.0**, not 0.1.2: an exported function gets new arguments, and
calibrating a winter crop with the new arguments gives substantially different
values. Under pre-1.0 semantic versioning, that is a minor bump.

Default `crop_type` is `"summer"` (as in the other functions), so a 0.1.1-style
call returns the same values as before (up to rounding). Users must pass
`crop_type = "winter"` and the winter rules to get the corrected winter
calibration. This keeps old calls stable but means a winter user who does not
read the NEWS gets the old behaviour (open question 2).

## 3. Impact on earlier results

Cross-version comparison: `validation/diagnostics/regression/` (each exported
model function run with 0.1.1 and 0.2.0 on identical inputs: bundled
historical weather and three DWD stations, 1991–2010, summer and winter
defaults). Log: `compare_v0.1.1_v0.2.0.log`.

| Exported function | Changes vs 0.1.1? | Example (old → new) |
| --- | --- | --- |
| `simulate_one_year()`, `run_simulation()` | No, for identical inputs (dates, statuses and accumulated TT identical) | DWD 00150 winter: 20 seasons, 14 mature → identical |
| `estimate_required_tt()` winter | No (bitwise identical) | 1604.1228947368 → 1604.1228947368 |
| `estimate_required_tt()` summer | Floating-point only, ≤ 1.2e-12 °C·day | 1267.8655250000 → 1267.8655250000 |
| `estimate_required_tt_from_observed()`, 0.1.1-style call | Floating-point only, ≤ 2.3e-12 °C·day; list has 2 extra elements | 1272.3712500000 → 1272.3712500000 |
| `estimate_required_tt_from_observed()` with `crop_type = "winter"` + rules | New behaviour (not possible in 0.1.1) | synthetic season 2048.82 → 1738.78; DWD wheat station medians ratio 0.779–0.876 (median 2378.7 → 1975.8) |
| `prepare_weather()`, `calc_daily_tt()`, `find_planting_date()`, `default_parameters()`, `compare_validation_metrics()`, `validate_weather_dwd()` | Not modified | – |

Consequences for results produced earlier with ThermoPheno:

- **Winter-crop requirements from `estimate_required_tt_from_observed()` in
  0.1.1** were raw sums (about 16% too high for the DWD wheat setup). Any
  winter-crop simulation driven by such a requirement tends to be late or
  forced at the harvest deadline. Such results should be recalibrated with
  `crop_type = "winter"` and the simulation's winter rules.
- **Requirements from `estimate_required_tt()`** (reference planting date +
  days to maturity) already applied the winter rules in 0.1.1 and are not
  affected. The Shiny app uses its own copy of this function (unchanged).
- **Summer crops**: values differ only by ~1e-12 °C·day. Because daily TT from
  DWD temperatures lies on a 0.05 °C·day grid, exact ties between a simulated
  total and the requirement occur, and rounding can then move a maturity date by
  one day. This happened once in the DWD validation (maize, station 10353, 2019:
  6 → 7 Sep, requirement 1099.275 → 1099.2750000000003).

## 4. Validation rerun

Script `validation/dwd_validate_thermopheno.R` with 0.2.0 installed. Outputs in
`pkgdown/assets/validation/`; 0.1.1 outputs archived unchanged in
`pkgdown/assets/validation/archive/v0.1.1/`. Comparison:
`validation/compare_validation_runs.R` and `validation/compare_common_cases.R`
(logs `comparison_v0.1.1_vs_new.log`, `comparison_common_cases.log`).

Inputs: all 48 DWD files have MD5 checksums identical to the 0.1.1 manifest. No
DWD archive changed. Same 75 attempted cases, same station matches, same
held-out station-years, identical simulated planting dates.

| Crop | Phase | n | Bias (d) | MAE (d) | RMSE (d) | R² |
| --- | --- | --- | --- | --- | --- | --- |
| Maize | Planting | 152 → 152 | 0.29 → 0.29 | 7.17 → 7.17 | 9.58 → 9.58 | 0.097 → 0.097 |
| Maize | Harvest | 152 → 152 | −7.34 → −7.34 | 15.33 → 15.32 | 19.58 → 19.58 | 0.251 → 0.252 |
| Sugar beet | Planting | 168 → 168 | −0.60 → −0.60 | 7.36 → 7.36 | 15.83 → 15.83 | 0.053 → 0.053 |
| Sugar beet | Harvest | 168 → 168 | −15.92 → −15.92 | 19.37 → 19.37 | 24.49 → 24.49 | 0.055 → 0.055 |
| Winter wheat | Planting | 217 → 217 | −8.76 → −8.76 | 10.39 → 10.39 | 13.49 → 13.49 | 0.002 → 0.002 |
| Winter wheat | Harvest | 209 → 217 | 13.33 → −1.79 | 14.52 → 9.45 | 17.18 → 11.12 | 0.096 → 0.269 |

Same 209 wheat station-years in both runs: bias 13.33 → −2.01, MAE 14.52 → 9.36,
RMSE 17.18 → 11.02, R² 0.096 → 0.259. Wheat harvest within 14 days: 55.0% → 84.3%.

| Winter wheat held-out status | 0.1.1 | 0.2.0 |
| --- | ---: | ---: |
| mature | 35 | 195 |
| forced_harvest_immature | 174 | 22 |
| failed_to_mature | 8 | 0 |
| insufficient_evaluation_weather | 2 | 2 |

Maize: 8 forced, 0 failed in both. Sugar beet: 0 forced, 0 failed in both.
Changed wheat harvests moved earlier by a median of 17 days (range 0 to −35).

Skill remains low for planting dates (R² ≤ 0.10) and sugar beet harvest. The
article states this plainly, keeps the 75-case reporting, the exclusion table,
the sugar beet 11908/2015 sensitivity note and the limitations, and adds
"What changed in 0.2.0" and "Comparison with 0.1.1" sections.

**Held-out data disclosure.** The hypothesis came from the 0.1.1 held-out
results (174 forced wheat harvests), which were already published. The
confirmation, the fix and the tests used only 1991–2010 and synthetic data,
and the deadline decision used only 1991–2010. I first saw 2011–2024 results
for 0.2.0 after the package code and the validation configuration were final.
After that I changed only reporting (comparison scripts, article, NEWS, this
report), not code or configuration that affects results.

## 5. Tests

New: `tests/testthat/test-tt-consistency.R` and `helper-synthetic-weather.R`.

- Round trip, winter and summer synthetic crops: calibrate from one season,
  simulate from the observed planting date, maturity within one day of the
  observed harvest (0.2.0: both exactly on the date).
- Raw mode: `accumulation = "raw"` equals the value produced by 0.1.1
  (2048.8223908266 winter, 1224.4860296517 summer) and the 0.1.1 formula, and
  fails the winter round trip (+19 days).
- Equality: shared accumulator ≡ simulation accumulation (`identical`), ≡ an
  independent re-implementation of the 0.1.1 simulation loop, ≡ calibration over
  the same interval; baseline calibration ≡ reference loop.
- Summer default ≡ raw up to 1e-12; winter simulation-rule calibration < raw.

Against 0.1.1: `validation/diagnostics/roundtrip/roundtrip_v0.1.1.log` shows the
winter round trip failing (maturity 13 Aug vs observed 25 Jul, +19 days) and the
summer one landing +1 day (rounding). `new_tests_on_v0.1.1.log`: the new test file
fails 6 tests on 0.1.1 (the new arguments do not exist there). On 0.2.0
(`roundtrip_v0.2.0.log`) both round trips are exact. All tests pass on 0.2.0
(69 expectations, 0 failures).

On real data, `validation/diagnostics/v0.2.0/diagnostic.log`: calibration with
`accumulation = "simulation"` equals the simulation's accumulation for all 1,047
calibration seasons (maximum difference 0).

## 6. R CMD check and CI

`rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))` on Windows,
R 4.5.2: **0 errors, 0 warnings, 2 NOTEs** (`validation/diagnostics/rcmdcheck_v0.2.0.log`).

- "CRAN incoming feasibility": New submission; tarball 11.9 MB. Pre-existing:
  the 0.1.1 tarball is 11,926,868 bytes vs 11,933,113 for 0.2.0. Mostly bundled
  example climate data (6.9 MB scenario CSV).
- "Non-standard top-level files": `RUN_CHECKS.R`, `SOURCE_PACKAGE_NOTES.md`,
  `scripts`, `validation`. Pre-existing. The new diagnostic files are excluded
  via `.Rbuildignore`.

Pandoc came from the local RStudio installation (`RSTUDIO_PANDOC`).
**GitHub CI has not run**: nothing has been pushed. pkgdown site built locally
without errors; all 27 local links in the validation article resolve.

## 7. Release details

Not done; waiting for your review of the changes on `main` (no PR, by your decision).

- Planned tag: `v0.2.0`, published as a GitHub release, which creates a new
  version under concept DOI **10.5281/zenodo.23619166**. `.zenodo.json` already
  says version 0.2.0 (MIT, same creator, description and related identifiers).
- Version DOI for 0.2.0: not yet assigned.
- After Zenodo assigns it: follow-up commit for README citation,
  `CITATION.cff` (version, DOI, `date-released`), `inst/CITATION` DOI map, and
  rebuild of the site. Until then `inst/CITATION` falls back to the concept DOI
  for 0.2.0; README and `CITATION.cff` still cite 0.1.1.
- Live article: https://mr-eini.github.io/ThermoPheno/articles/dwd-validation.html
  (updates after merge). Archived 0.1.1 results (after merge):
  https://github.com/MR-Eini/ThermoPheno/tree/main/pkgdown/assets/validation/archive/v0.1.1
  and on the site under `validation/archive/v0.1.1/`. 0.1.1 article source:
  https://github.com/MR-Eini/ThermoPheno/blob/v0.1.1/vignettes/articles/dwd-validation.Rmd

## 8. Open questions for the author

1. **Winter wheat harvest deadline (15 Aug).** 11% of calibration-period
   harvests at the selected stations were later. 22 held-out harvests are still
   forced. I did not change it, and since I have now seen 2011–2024 results, any
   change should rest on a documented agronomic source chosen independently of
   them.
2. **Default `crop_type = "summer"` in `estimate_required_tt_from_observed()`.**
   It keeps old calls stable, but a winter user who omits `crop_type` gets the
   raw-equivalent sum. Alternatives: warn when the call looks like a winter
   crop, or require `crop_type` (breaking).
3. **Maturity ties.** The 0.05 °C·day grid makes exact ties common, so maturity
   on tie days depends on rounding. A small tolerance (e.g. `cum >= req - 1e-9`)
   would make it robust but changes simulation output in rare cases.
4. **Shiny app** (`inst/app/app.R`) keeps its own copies of the model functions.
   It does not use observed calibration, so it is unaffected, but it is a second
   source of truth. Refactoring it to call the package functions is a separate
   task.
5. **Vernalization temperature window** is not in `crop_config.csv`; the defaults
   (0–10 °C) are used on both sides. You may want it recorded explicitly.
6. **Branch.** You decided to keep the changes on `main` without a separate
   branch or pull request. Review of the working-tree changes therefore stands
   in for PR approval before the `v0.2.0` release. CI runs once the commit is
   pushed to `main`.

## Environment notes

- ThermoPheno 0.2.0 (this working tree) was installed into the user R library
  (`C:/Users/p101727/AppData/Local/R/win-library/4.5`), replacing 0.1.0 there.
  0.1.1 was installed in a temporary library only for the baseline runs.
- DWD inputs are cached in `validation/cache/` (git-ignored).
