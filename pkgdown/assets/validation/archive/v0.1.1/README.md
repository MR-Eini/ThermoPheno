# Archived DWD validation outputs: ThermoPheno 0.1.1

These files are the unchanged outputs of the DWD validation published with
ThermoPheno 0.1.1 (git tag `v0.1.1`). They are kept so that the 0.1.1 and
later results can be compared. The current outputs are one level up, in
`../../`.

In 0.1.1, the winter wheat thermal-time requirement was calibrated as a raw
triangular thermal-time sum over the observed planting–harvest interval,
without the vernalization and dormancy rules that the simulation applies.
ThermoPheno 0.2.0 calibrates with the simulation rules. See the validation
article for the comparison.

To reproduce these files, install the `v0.1.1` tag and run
`validation/dwd_validate_thermopheno.R` from that tag:

```sh
Rscript validation/dwd_validate_thermopheno.R
```

Article source for 0.1.1:
<https://github.com/MR-Eini/ThermoPheno/blob/v0.1.1/vignettes/articles/dwd-validation.Rmd>

Source: Deutscher Wetterdienst (DWD), Climate Data Center.
[DWD CDC terms of use](https://opendata.dwd.de/climate_environment/CDC/Terms_of_use.pdf) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
