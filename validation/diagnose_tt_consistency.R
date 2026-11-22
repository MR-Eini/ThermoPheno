# Calibration/simulation thermal-time consistency diagnostic.
#
# Uses CALIBRATION harvest years (1991-2010) only. No 2011-2024 observation is
# read into the comparison tables. The DWD input discovery, download and parsing
# code is taken unchanged from validation/dwd_validate_thermopheno.R (everything
# before the per-crop loop is evaluated from that file), and the station
# selection below copies the rule used there (25 stations per crop with the most
# planting-harvest calendars in 1991-2010, ties broken by station ID).
#
# For every calibration station-year with an observed planting and harvest date:
#   (a) raw_tt: estimate_required_tt_from_observed(), i.e. the triangular sum over
#       the observed interval exactly as the calibration does it;
#   (b) sim_tt: the thermal time simulate_one_year() accumulates over the same
#       interval when planting is forced to the observed date (planting window
#       collapsed to that day, temperature gates disabled, requirement = Inf, no
#       forced harvest, season end = observed harvest date). This calls the
#       simulation code path itself, so (b) is what the model really integrates.
# and, for each station:
#   (c) reach_tt: thermal time the simulation accumulates from its own simulated
#       planting date to the configured harvest deadline, per calibration year,
#       compared with the calibrated requirement.
#
# Run from repository root with the package version to diagnose installed:
#   Rscript validation/diagnose_tt_consistency.R [output_dir]

args <- commandArgs(trailingOnly = TRUE)
diag_dir <- if (length(args) >= 1) args[[1]] else "validation/diagnostics"
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Evaluate the shared set-up section of the validation script ------------
val_exprs <- parse("validation/dwd_validate_thermopheno.R", keep.source = FALSE)
stop_at <- which(vapply(val_exprs, function(e) {
  identical(deparse(e)[1], "all_pairs <- list()")
}, logical(1)))
if (length(stop_at) != 1L) stop("Could not locate the per-crop loop in the validation script.")
for (e in val_exprs[seq_len(stop_at - 1L)]) eval(e, envir = globalenv())

message("ThermoPheno version under diagnosis: ", packageVersion("ThermoPheno"))
message("Calibration years only: ", min(calib_years), "-", max(calib_years))

season_rows <- list()
reach_rows <- list()
harvest_obs_rows <- list()

for (ii in seq_len(nrow(cfg))) {
  cc <- cfg[ii, ]
  is_winter_crop <- identical(as.character(cc$crop_type[[1]]), "winter")
  message("Diagnosing crop: ", cc$crop_name)

  # ---- Observation pairing: copied from the validation script -------------
  chosen <- pheno_files %>% filter(file == cc$phenology_file)
  if (nrow(chosen) != 1L) stop("Could not locate DWD input for ", cc$crop_name)
  pheno_path <- download_cached(paste0(chosen$base_url, chosen$file))
  pheno_raw <- read_dwd_table(pheno_path)
  pheno <- join_phase_names(pheno_raw, phase_def) %>%
    filter(object_id_tmp == as.character(cc$object_id))
  pheno <- pheno %>% filter(is.na(qb_tmp) | !(qb_tmp %in% c("5", "8")))
  planting_obs <- pheno %>% filter(phase_id_tmp == "10", !is.na(date_obs))
  harvest_obs  <- pheno %>% filter(phase_id_tmp == "24", !is.na(date_obs))
  planting_obs <- planting_obs %>% mutate(
    crop_year = dplyr::if_else(
      .env$is_winter_crop & lubridate::month(date_obs) >= 8L,
      lubridate::year(date_obs) + 1L,
      lubridate::year(date_obs)
    )
  )
  harvest_obs <- harvest_obs %>% mutate(crop_year = lubridate::year(date_obs))
  median_date <- function(x) {
    x <- as.Date(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(as.Date(NA))
    as.Date(stats::median(as.numeric(x)), origin = "1970-01-01")
  }
  planting_year <- planting_obs %>%
    group_by(station_id_pheno, crop_year) %>%
    summarise(observed_planting_date = median_date(date_obs), .groups = "drop")
  harvest_year <- harvest_obs %>%
    group_by(station_id_pheno, crop_year) %>%
    summarise(observed_harvest_date = median_date(date_obs), .groups = "drop")
  # Diagnostic restriction: calibration years only.
  obs <- inner_join(planting_year, harvest_year,
                    by = c("station_id_pheno", "crop_year"),
                    relationship = "one-to-one") %>%
    filter(crop_year %in% calib_years,
           observed_harvest_date >= observed_planting_date)

  selected_sites <- obs %>%
    count(station_id_pheno, name = "calibration_calendar_years") %>%
    filter(calibration_calendar_years >= MIN_OBSERVED_TT_YEARS) %>%
    inner_join(pheno_sites_all, by = "station_id_pheno", relationship = "many-to-one") %>%
    arrange(desc(calibration_calendar_years), station_id_pheno) %>%
    head(25)

  harvest_obs_rows[[cc$crop_name]] <- obs %>%
    filter(station_id_pheno %in% selected_sites$station_id_pheno) %>%
    transmute(crop_name = cc$crop_name, station_id_pheno, crop_year,
              observed_planting_date, observed_harvest_date)

  for (ss in seq_len(nrow(selected_sites))) {
    site <- selected_sites[ss, ]
    d <- haversine_km(site$pheno_lon, site$pheno_lat, clim_meta$lon, clim_meta$lat)
    nearest <- clim_meta[which.min(d), , drop = FALSE]
    obs_site <- obs %>% filter(station_id_pheno == site$station_id_pheno)
    weather <- tryCatch(read_climate_station(nearest), error = function(e) NULL)
    if (is.null(weather)) next

    # (a) raw calibration sums, exactly as the 0.1.1 validation calibrated.
    cal_args <- list(
      weather = weather, observed_calendar = obs_site,
      calibration_years = calib_years,
      planting_col = "observed_planting_date",
      harvest_col = "observed_harvest_date", year_col = "crop_year",
      t_base = cc$t_base, t_opt = cc$t_opt, t_max_cut = cc$t_max_cut,
      tt_mode = cc$tt_mode, summary_fun = OBSERVED_TT_SUMMARY_FUN,
      min_weather_coverage = MIN_WEATHER_COVERAGE
    )
    has_accumulation <- "accumulation" %in% names(formals(ThermoPheno::estimate_required_tt_from_observed))
    req <- tryCatch(
      do.call(ThermoPheno::estimate_required_tt_from_observed,
              c(cal_args, if (has_accumulation) list(accumulation = "raw"))),
      error = function(e) NULL
    )
    # (a') with ThermoPheno >= 0.2.0: the calibration with simulation rules.
    req_sim <- if (has_accumulation) tryCatch(
      do.call(ThermoPheno::estimate_required_tt_from_observed, c(cal_args, list(
        accumulation = "simulation", crop_type = cc$crop_type,
        winter_dormancy_temp = cc$winter_dormancy_temp,
        vernalization_required = isTRUE(cc$vernalization_required),
        vernalization_days_required = cc$vernalization_days_required,
        spring_regrowth_temp = cc$spring_regrowth_temp
      ))),
      error = function(e) NULL
    ) else NULL
    if (is.null(req)) next
    yr_tt <- req$yearly_required_tt
    yr_tt <- yr_tt[is.finite(yr_tt$observed_required_tt), , drop = FALSE]
    if (!nrow(yr_tt)) next
    calibrated_ok <- req$n_valid_years >= MIN_OBSERVED_TT_YEARS &&
      is.finite(req$required_tt) && req$required_tt > 0

    # (b) simulation accumulation over the same observed interval.
    sim_tt <- vapply(seq_len(nrow(yr_tt)), function(k) {
      p <- yr_tt$observed_planting_date[k]
      h <- yr_tt$observed_harvest_date[k]
      expected_plant_year <- yr_tt$crop_year[k] - as.integer(is_winter_crop)
      if (as.integer(format(p, "%Y")) != expected_plant_year ||
          as.integer(format(h, "%Y")) != yr_tt$crop_year[k]) return(NA_real_)
      r <- ThermoPheno::simulate_one_year(
        weather = weather, sim_year = yr_tt$crop_year[k], required_tt = Inf,
        earliest_planting_mmdd = format(p, "%m-%d"),
        latest_planting_mmdd = format(p, "%m-%d"),
        latest_harvest_mmdd = format(h, "%m-%d"),
        t_base = cc$t_base, t_opt = cc$t_opt, t_max_cut = cc$t_max_cut,
        tt_mode = cc$tt_mode, crop_type = cc$crop_type,
        min_mean_temp_plant = -999, forced_harvest_allowed = FALSE,
        winter_dormancy_temp = cc$winter_dormancy_temp,
        vernalization_required = isTRUE(cc$vernalization_required),
        vernalization_days_required = cc$vernalization_days_required,
        spring_regrowth_temp = cc$spring_regrowth_temp,
        winter_plant_temp_min = -999, winter_plant_temp_max = 999
      )
      if (!identical(as.Date(r$planting_date), p)) return(NA_real_)
      r$accumulated_tt
    }, numeric(1))

    season_rows[[paste(cc$crop_name, site$station_id_pheno)]] <- data.frame(
      crop_name = cc$crop_name, crop_type = cc$crop_type,
      station_id_pheno = site$station_id_pheno,
      climate_station_id = nearest$station_id,
      crop_year = yr_tt$crop_year,
      observed_planting_date = yr_tt$observed_planting_date,
      observed_harvest_date = yr_tt$observed_harvest_date,
      raw_tt = yr_tt$observed_required_tt,
      sim_tt = sim_tt,
      ratio_sim_to_raw = sim_tt / yr_tt$observed_required_tt,
      station_calibrated = calibrated_ok,
      calibrated_required_tt_raw = req$required_tt,
      calibration_sim_tt = if (is.null(req_sim)) NA_real_ else {
        ys <- req_sim$yearly_required_tt
        ys$observed_required_tt[match(yr_tt$crop_year, ys$crop_year)]
      }
    )

    # (c) thermal time reachable before the harvest deadline, simulation
    #     planting rules as configured, calibration years only.
    if (calibrated_ok) {
      cal_years_wx <- intersect(calib_years, unique(weather$year))
      if (is_winter_crop) cal_years_wx <- cal_years_wx[(cal_years_wx - 1L) %in% weather$year]
      reach <- lapply(cal_years_wx, function(yr) {
        r <- ThermoPheno::simulate_one_year(
          weather = weather, sim_year = yr, required_tt = Inf,
          earliest_planting_mmdd = cc$earliest_planting_mmdd,
          latest_planting_mmdd = cc$latest_planting_mmdd,
          latest_harvest_mmdd = cc$latest_harvest_mmdd,
          t_base = cc$t_base, t_opt = cc$t_opt, t_max_cut = cc$t_max_cut,
          tt_mode = cc$tt_mode, crop_type = cc$crop_type,
          min_mean_temp_plant = cc$min_mean_temp_plant,
          forced_harvest_allowed = FALSE,
          winter_dormancy_temp = cc$winter_dormancy_temp,
          vernalization_required = isTRUE(cc$vernalization_required),
          vernalization_days_required = cc$vernalization_days_required,
          spring_regrowth_temp = cc$spring_regrowth_temp,
          winter_plant_temp_min = cc$winter_plant_temp_min,
          winter_plant_temp_max = cc$winter_plant_temp_max
        )
        data.frame(crop_year = yr, simulated_planting_date = r$planting_date,
                   status = r$status, reach_tt = r$accumulated_tt)
      })
      reach <- do.call(rbind, reach)
      reach$crop_name <- cc$crop_name
      reach$station_id_pheno <- site$station_id_pheno
      reach$calibrated_required_tt_raw <- req$required_tt
      reach$calibrated_required_tt_sim <- stats::median(sim_tt, na.rm = TRUE)
      reach_rows[[paste(cc$crop_name, site$station_id_pheno)]] <- reach
    }
  }
}

seasons <- bind_rows(season_rows)
reach <- bind_rows(reach_rows)
harvests <- bind_rows(harvest_obs_rows)
stopifnot(all(seasons$crop_year %in% calib_years), all(reach$crop_year %in% calib_years))

q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
ratio_summary <- function(df) {
  df %>% summarise(
    n_seasons = sum(is.finite(ratio_sim_to_raw)),
    n_excluded = sum(!is.finite(ratio_sim_to_raw)),
    ratio_min = suppressWarnings(min(ratio_sim_to_raw, na.rm = TRUE)),
    ratio_p10 = q(ratio_sim_to_raw, 0.10),
    ratio_median = q(ratio_sim_to_raw, 0.50),
    ratio_mean = mean(ratio_sim_to_raw, na.rm = TRUE),
    ratio_p90 = q(ratio_sim_to_raw, 0.90),
    ratio_max = suppressWarnings(max(ratio_sim_to_raw, na.rm = TRUE)),
    share_ratio_below_0_99 = mean(ratio_sim_to_raw < 0.99, na.rm = TRUE),
    median_raw_minus_sim_tt = q(raw_tt - sim_tt, 0.5),
    .groups = "drop"
  )
}
by_station <- seasons %>% group_by(crop_name, station_id_pheno, climate_station_id) %>% ratio_summary()
by_crop <- seasons %>% group_by(crop_name) %>% ratio_summary()

reach_summary <- reach %>%
  group_by(crop_name) %>%
  summarise(
    n_station_years = n(),
    n_stations = n_distinct(station_id_pheno),
    n_not_planted = sum(status == "not_planted"),
    share_reach_below_raw_requirement = mean(reach_tt < calibrated_required_tt_raw, na.rm = TRUE),
    share_reach_below_sim_requirement = mean(reach_tt < calibrated_required_tt_sim, na.rm = TRUE),
    median_reach_tt = q(reach_tt, 0.5),
    median_raw_requirement = q(calibrated_required_tt_raw, 0.5),
    median_sim_requirement = q(calibrated_required_tt_sim, 0.5),
    .groups = "drop"
  )

deadline_check <- harvests %>%
  left_join(cfg %>% select(crop_name, latest_harvest_mmdd), by = "crop_name") %>%
  mutate(deadline = as.Date(paste(crop_year, latest_harvest_mmdd, sep = "-")),
         harvest_doy = lubridate::yday(observed_harvest_date)) %>%
  group_by(crop_name, latest_harvest_mmdd) %>%
  summarise(
    n_calibration_harvests = n(),
    share_observed_after_deadline = mean(observed_harvest_date > deadline),
    harvest_p50 = format(as.Date(q(harvest_doy, 0.5) - 1, origin = "2001-01-01"), "%m-%d"),
    harvest_p90 = format(as.Date(q(harvest_doy, 0.9) - 1, origin = "2001-01-01"), "%m-%d"),
    harvest_p95 = format(as.Date(q(harvest_doy, 0.95) - 1, origin = "2001-01-01"), "%m-%d"),
    harvest_max = format(as.Date(max(harvest_doy) - 1, origin = "2001-01-01"), "%m-%d"),
    .groups = "drop"
  )

readr::write_csv(seasons, file.path(diag_dir, "tt_consistency_by_season.csv"))
readr::write_csv(by_station, file.path(diag_dir, "tt_consistency_by_station.csv"))
readr::write_csv(by_crop, file.path(diag_dir, "tt_consistency_by_crop.csv"))
readr::write_csv(reach, file.path(diag_dir, "tt_reachable_before_deadline_by_year.csv"))
readr::write_csv(reach_summary, file.path(diag_dir, "tt_reachable_before_deadline_summary.csv"))
readr::write_csv(deadline_check, file.path(diag_dir, "calibration_harvest_vs_deadline.csv"))

options(width = 200)
cat("ThermoPheno version:", as.character(packageVersion("ThermoPheno")), "\n")
cat("Calibration years:", min(calib_years), "-", max(calib_years), "\n\n")
cat("(b)/(a) by crop\n"); print(as.data.frame(by_crop))
if (any(is.finite(seasons$calibration_sim_tt))) {
  cat("\nCalibration with accumulation = 'simulation' vs simulation accumulation (b)\n")
  print(as.data.frame(seasons %>% group_by(crop_name) %>% summarise(
    n = sum(is.finite(calibration_sim_tt) & is.finite(sim_tt)),
    max_abs_difference = max(abs(calibration_sim_tt - sim_tt), na.rm = TRUE),
    .groups = "drop")))
}
cat("\nReachable thermal time before the harvest deadline\n"); print(as.data.frame(reach_summary))
cat("\nObserved calibration-period harvest dates vs configured deadline\n"); print(as.data.frame(deadline_check))
