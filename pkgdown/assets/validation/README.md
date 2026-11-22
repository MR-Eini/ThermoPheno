# DWD phenology validation

Install ThermoPheno 0.2.0 from this repository. The script loads the installed
package and uses `readr`, `dplyr`, `lubridate`, `ggplot2` and `tidyr`.
From the repository root, run:

```sh
Rscript validation/dwd_validate_thermopheno.R
```

Inputs are cached in `validation/cache/` and outputs go to
`validation/results/`. The input manifest records URLs, cache file dates and
MD5 checksums. Use an empty cache to download fresh copies. DWD can update its
archives. Compare the manifest with the one shipped with the website when
reproducing the published tables and figures. Site assets are kept in
`pkgdown/assets/validation/`.

Source: Deutscher Wetterdienst (DWD), Climate Data Center. The workflow filters,
pairs and aggregates the observations and compares the result with simulations.
[DWD CDC terms of use](https://opendata.dwd.de/climate_environment/CDC/Terms_of_use.pdf) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
