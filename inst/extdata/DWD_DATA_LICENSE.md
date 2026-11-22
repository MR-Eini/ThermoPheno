# DWD observed weather example

Source: Deutscher Wetterdienst (DWD), Climate Data Center. Daily station
observations for Germany, historical KL product, station 04709
(Simmern-Wahlbach). [Source archive](https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_04709_19790401_20251231_hist.zip)

The source archive was accessed on 22 November 2026. Its SHA-256 is
`324a5545c329c9fb08944ad5a8e01064a55c285e4a3f124c2c12d0a64e045fb0`.
The CSV selects 1 January through 31 December 2005, reformats the date,
renames TNK and TXK to `tmin` and `tmax` (°C), and preserves the station ID.
The temperatures are unchanged. See `validation/extract_dwd_example.py` in
the repository for the extraction.

The DWD data are licensed under Creative Commons Attribution 4.0 International.
The ThermoPheno software has its own MIT license. [DWD CDC terms of use](https://opendata.dwd.de/climate_environment/CDC/Terms_of_use.pdf) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
