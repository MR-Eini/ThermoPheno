# Side-by-side comparison of two DWD validation runs (metrics, harvest outcomes,
# exclusions and input checksums).
#   Rscript validation/compare_validation_runs.R <old_dir> <new_dir> <out_dir>
# Defaults: archived 0.1.1 outputs vs validation/results.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })
a <- commandArgs(trailingOnly = TRUE)
old_dir <- if (length(a) >= 1) a[[1]] else "pkgdown/assets/validation/archive/v0.1.1"
new_dir <- if (length(a) >= 2) a[[2]] else "validation/results"
out_dir <- if (length(a) >= 3) a[[3]] else "validation/results"
rd <- function(d, f) read_csv(file.path(d, f), show_col_types = FALSE,
                              col_types = cols(station_id_pheno = col_character(),
                                               climate_station_id = col_character(),
                                               .default = col_guess()))

m_old <- rd(old_dir, "validation_metrics.csv"); m_new <- rd(new_dir, "validation_metrics.csv")
metrics_cmp <- full_join(
  m_old %>% select(crop_name, phase, n, bias, mae, rmse, r2),
  m_new %>% select(crop_name, phase, n, bias, mae, rmse, r2),
  by = c("crop_name", "phase"), suffix = c("_v011", "_new")
)

p_old <- rd(old_dir, "validation_pairs.csv"); p_new <- rd(new_dir, "validation_pairs.csv")
status_cmp <- full_join(
  p_old %>% count(crop_name, status, name = "n_v011"),
  p_new %>% count(crop_name, status, name = "n_new"),
  by = c("crop_name", "status")
) %>% mutate(across(starts_with("n_"), ~ coalesce(.x, 0L))) %>% arrange(crop_name, status)

same_cases <- nrow(p_old) == nrow(p_new) &&
  identical(paste(p_old$crop_name, p_old$station_id_pheno, p_old$year),
            paste(p_new$crop_name, p_new$station_id_pheno, p_new$year))
planting_identical <- same_cases &&
  identical(as.character(p_old$simulated_planting_date), as.character(p_new$simulated_planting_date))

x_old <- rd(old_dir, "nearest_climate_station_matches.csv")
x_new <- rd(new_dir, "nearest_climate_station_matches.csv")
exclusions <- full_join(
  x_old %>% count(crop_name, evaluation_status, name = "n_v011"),
  x_new %>% count(crop_name, evaluation_status, name = "n_new"),
  by = c("crop_name", "evaluation_status")
) %>% mutate(across(starts_with("n_"), ~ coalesce(.x, 0L)))
matches_identical <- identical(
  x_old %>% select(crop_name, station_id_pheno, climate_station_id, distance_km),
  x_new %>% select(crop_name, station_id_pheno, climate_station_id, distance_km))

i_old <- rd(old_dir, "dwd_input_manifest.csv"); i_new <- rd(new_dir, "dwd_input_manifest.csv")
inputs <- full_join(i_old %>% select(url, md5_v011 = md5), i_new %>% select(url, md5_new = md5), by = "url") %>%
  mutate(check = case_when(is.na(md5_v011) ~ "new_input",
                           is.na(md5_new) ~ "not_used_now",
                           md5_v011 == md5_new ~ "identical",
                           TRUE ~ "changed"))

r_old <- rd(old_dir, "required_tt_calibration_by_station.csv")
r_new <- rd(new_dir, "required_tt_calibration_by_station.csv")
req_cmp <- full_join(
  r_old %>% distinct(crop_name, station_id_pheno, required_tt_v011 = required_tt_final),
  r_new %>% distinct(crop_name, station_id_pheno, required_tt_new = required_tt_final),
  by = c("crop_name", "station_id_pheno")
) %>% mutate(ratio_new_to_v011 = required_tt_new / required_tt_v011)
req_summary <- req_cmp %>% group_by(crop_name) %>% summarise(
  n_stations = n(),
  median_required_tt_v011 = median(required_tt_v011),
  median_required_tt_new = median(required_tt_new),
  min_ratio = min(ratio_new_to_v011), median_ratio = median(ratio_new_to_v011),
  max_ratio = max(ratio_new_to_v011), .groups = "drop")

write_csv(metrics_cmp, file.path(out_dir, "comparison_metrics_v0.1.1_vs_new.csv"))
write_csv(status_cmp, file.path(out_dir, "comparison_status_counts_v0.1.1_vs_new.csv"))
write_csv(exclusions, file.path(out_dir, "comparison_exclusions_v0.1.1_vs_new.csv"))
write_csv(inputs, file.path(out_dir, "comparison_input_checksums_v0.1.1_vs_new.csv"))
write_csv(req_cmp, file.path(out_dir, "comparison_required_tt_v0.1.1_vs_new.csv"))

options(width = 200, pillar.sigfig = 4)
cat("Metrics\n"); print(as.data.frame(metrics_cmp), digits = 4)
cat("\nHeld-out outcome status counts\n"); print(as.data.frame(status_cmp))
cat("\nSame held-out station-year cases in both runs:", same_cases, "\n")
cat("Simulated planting dates identical in both runs:", planting_identical, "\n")
cat("Station matching identical:", matches_identical, "\n")
cat("\nAttempted-case outcomes\n"); print(as.data.frame(exclusions))
cat("\nInput checksums\n"); print(table(inputs$check))
if (any(inputs$check != "identical")) print(as.data.frame(filter(inputs, check != "identical")))
cat("\nCalibrated requirement by crop\n"); print(as.data.frame(req_summary), digits = 5)
