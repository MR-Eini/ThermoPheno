packages <- c(
  "shiny",
  "dplyr",
  "lubridate",
  "ggplot2",
  "DT",
  "readr",
  "ggridges",
  "RColorBrewer",
  "tidyr",
  "tibble"
)

missing_pkgs <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
invisible(lapply(packages, library, character.only = TRUE))

# ============================================================
# PhenoFlex - Historical + Climate Change app
# Uses historical data to estimate required thermal time
# Then runs historical and climate simulations with same crop parameters
# No frost logic
# No user-defined planting/harvest windows in the UI
# Uses only maximum growing days
# ============================================================

safe_as_date <- function(x) as.Date(x)

options(shiny.maxRequestSize = 200 * 1024^2)

clamp_year_day <- function(year, mmdd) {
  as.Date(paste0(year, "-", mmdd))
}

prepare_weather <- function(df) {
  names(df) <- tolower(names(df))
  required <- c("date", "tmin", "tmax")
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(paste("Missing required weather columns:", paste(miss, collapse = ", ")))
  
  df %>%
    mutate(
      date = safe_as_date(date),
      year = year(date),
      doy = yday(date),
      tmean = (tmin + tmax) / 2
    ) %>%
    arrange(date)
}

calc_daily_tt <- function(tmin, tmax, t_base, t_opt = NA, t_max_cut = NA,
                          mode = c("simple", "capped", "triangular")) {
  mode <- match.arg(mode)
  tmean <- (tmin + tmax) / 2
  
  if (mode == "simple") {
    tt <- pmax(tmean - t_base, 0)
  } else if (mode == "capped") {
    if (is.na(t_opt)) t_opt <- 999
    tt <- pmax(pmin(tmean, t_opt) - t_base, 0)
  } else {
    if (is.na(t_opt) || is.na(t_max_cut)) stop("For triangular mode, both t_opt and t_max_cut must be provided.")
    tt <- numeric(length(tmean))
    idx1 <- tmean <= t_base
    idx2 <- tmean > t_base & tmean <= t_opt
    idx3 <- tmean > t_opt & tmean < t_max_cut
    idx4 <- tmean >= t_max_cut
    tt[idx1] <- 0
    tt[idx2] <- tmean[idx2] - t_base
    peak <- t_opt - t_base
    tt[idx3] <- peak * (t_max_cut - tmean[idx3]) / (t_max_cut - t_opt)
    tt[idx4] <- 0
    tt[tt < 0] <- 0
  }
  
  tt
}

default_parameters <- function(crop_type = c("summer", "winter")) {
  crop_type <- match.arg(crop_type)
  if (crop_type == "summer") {
    list(
      crop_name = "Maize",
      days_to_maturity = 140,
      t_base = 8,
      t_opt = 25,
      t_max_cut = 35,
      baseline_planting_mmdd = "04-15",
      max_growing_days = 220,
      min_mean_temp_plant = 8,
      winter_plant_temp_min = 5,
      winter_plant_temp_max = 15
    )
  } else {
    list(
      crop_name = "Winter wheat",
      days_to_maturity = 300,
      t_base = 0,
      t_opt = 18,
      t_max_cut = 30,
      baseline_planting_mmdd = "10-01",
      max_growing_days = 320,
      min_mean_temp_plant = 8,
      winter_plant_temp_min = 5,
      winter_plant_temp_max = 15
    )
  }
}

get_internal_windows <- function(crop_type = c("summer", "winter")) {
  crop_type <- match.arg(crop_type)
  if (crop_type == "summer") {
    list(
      earliest_planting_mmdd = "03-01",
      latest_planting_mmdd   = "06-30"
    )
  } else {
    list(
      earliest_planting_mmdd = "09-01",
      latest_planting_mmdd   = "11-30"
    )
  }
}

estimate_required_tt <- function(weather,
                                 baseline_years,
                                 planting_mmdd,
                                 days_to_maturity,
                                 t_base,
                                 t_opt = NA,
                                 t_max_cut = NA,
                                 tt_mode = "simple",
                                 crop_type = c("summer", "winter"),
                                 winter_dormancy_temp = 0,
                                 vernalization_required = FALSE,
                                 vernalization_temp_min = 0,
                                 vernalization_temp_max = 10,
                                 vernalization_days_required = 30,
                                 spring_regrowth_temp = 5) {
  crop_type <- match.arg(crop_type)
  yrs <- sort(unique(weather$year))
  yrs <- yrs[yrs %in% baseline_years]
  if (length(yrs) == 0) stop("No overlap between baseline years and weather years.")
  
  if (crop_type == "summer") {
    yearly <- lapply(yrs, function(yr) {
      plant_date <- clamp_year_day(yr, planting_mmdd)
      end_date <- plant_date + days(days_to_maturity - 1)
      sub <- weather %>% filter(date >= plant_date, date <= end_date)
      if (nrow(sub) < days_to_maturity) return(data.frame(year = yr, required_tt = NA_real_))
      tt <- calc_daily_tt(sub$tmin, sub$tmax, t_base, t_opt, t_max_cut, mode = tt_mode)
      data.frame(year = yr, required_tt = sum(tt, na.rm = TRUE))
    }) %>% bind_rows()
  } else {
    yearly <- lapply(yrs, function(yr) {
      plant_date <- clamp_year_day(yr, planting_mmdd)
      end_date <- plant_date + days(days_to_maturity - 1)
      sub <- weather %>% filter(date >= plant_date, date <= end_date) %>% arrange(date)
      if (nrow(sub) < days_to_maturity) return(data.frame(year = yr, required_tt = NA_real_))
      
      vernal_days <- 0
      vern_sat <- ifelse(vernalization_required, FALSE, TRUE)
      regrowth_started <- FALSE
      cum_tt <- 0
      
      for (i in seq_len(nrow(sub))) {
        tmean_i <- sub$tmean[i]
        if (vernalization_required && !vern_sat) {
          if (tmean_i >= vernalization_temp_min && tmean_i <= vernalization_temp_max) vernal_days <- vernal_days + 1
          if (vernal_days >= vernalization_days_required) vern_sat <- TRUE
        }
        
        if (vern_sat && !regrowth_started && tmean_i >= spring_regrowth_temp) regrowth_started <- TRUE
        
        if (!vern_sat) {
          tt_i <- 0
        } else {
          tt_i <- calc_daily_tt(sub$tmin[i], sub$tmax[i], t_base, t_opt, t_max_cut, mode = tt_mode)
          if (tmean_i <= winter_dormancy_temp) tt_i <- 0
        }
        
        cum_tt <- cum_tt + tt_i
      }
      data.frame(year = yr, required_tt = cum_tt)
    }) %>% bind_rows()
  }
  
  valid <- yearly %>% filter(is.finite(required_tt))
  if (nrow(valid) == 0) stop("Could not estimate required thermal time from baseline period.")
  
  list(
    yearly_required_tt = yearly,
    required_tt = mean(valid$required_tt, na.rm = TRUE)
  )
}

find_planting_date <- function(weather_window,
                               earliest_planting_date,
                               latest_planting_date,
                               crop_type = c("summer", "winter"),
                               min_mean_temp_plant = -999,
                               winter_plant_temp_min = 5,
                               winter_plant_temp_max = 15) {
  crop_type <- match.arg(crop_type)
  
  candidates <- weather_window %>%
    filter(date >= earliest_planting_date, date <= latest_planting_date) %>%
    arrange(date)
  
  if (nrow(candidates) == 0) return(list(found = FALSE, planting_date = as.Date(NA)))
  
  for (i in seq_len(nrow(candidates))) {
    this_date <- candidates$date[i]
    tmean_today <- candidates$tmean[i]
    
    if (crop_type == "summer") {
      condition <- tmean_today >= min_mean_temp_plant
    } else {
      condition <- tmean_today >= winter_plant_temp_min && tmean_today <= winter_plant_temp_max
    }
    
    if (condition) return(list(found = TRUE, planting_date = this_date))
  }
  
  list(found = FALSE, planting_date = as.Date(NA))
}

simulate_one_year <- function(weather,
                              sim_year,
                              required_tt,
                              max_growing_days,
                              t_base,
                              t_opt = NA,
                              t_max_cut = NA,
                              tt_mode = "simple",
                              crop_type = c("summer", "winter"),
                              min_mean_temp_plant = 0,
                              forced_harvest_allowed = TRUE,
                              min_fraction_tt_for_forced_harvest = 0,
                              winter_dormancy_temp = 0,
                              vernalization_required = FALSE,
                              vernalization_temp_min = 0,
                              vernalization_temp_max = 10,
                              vernalization_days_required = 30,
                              spring_regrowth_temp = 5,
                              winter_plant_temp_min = 5,
                              winter_plant_temp_max = 15,
                              crop_name = NA_character_) {
  crop_type <- match.arg(crop_type)
  yr <- sim_year
  windows <- get_internal_windows(crop_type)
  
  if (crop_type == "summer") {
    weather_window <- weather %>% filter(year == yr)
    earliest_planting_date <- clamp_year_day(yr, windows$earliest_planting_mmdd)
    latest_planting_date   <- clamp_year_day(yr, windows$latest_planting_mmdd)
  } else {
    weather_window <- weather %>%
      filter(date >= clamp_year_day(yr - 1, windows$earliest_planting_mmdd),
             date <= clamp_year_day(yr, "12-31"))
    earliest_planting_date <- clamp_year_day(yr - 1, windows$earliest_planting_mmdd)
    latest_planting_date   <- clamp_year_day(yr - 1, windows$latest_planting_mmdd)
  }
  
  plant_result <- find_planting_date(
    weather_window = weather_window,
    earliest_planting_date = earliest_planting_date,
    latest_planting_date = latest_planting_date,
    crop_type = crop_type,
    min_mean_temp_plant = min_mean_temp_plant,
    winter_plant_temp_min = winter_plant_temp_min,
    winter_plant_temp_max = winter_plant_temp_max
  )
  
  if (!plant_result$found) {
    return(data.frame(
      crop_name = crop_name, crop_type = crop_type, row_role = "harvest_year", year = yr,
      planting_date = as.Date(NA), maturity_date = as.Date(NA), harvest_date = as.Date(NA),
      season_length_days = NA_integer_, accumulated_tt = 0, required_tt = required_tt,
      maturity_fraction = 0, status = "not_planted", forced_harvest = FALSE,
      regrowth_started = FALSE,
      vernalization_satisfied = ifelse(crop_type == "winter" && vernalization_required, FALSE, NA),
      vernalization_days = NA_real_
    ))
  }
  
  planting_date <- plant_result$planting_date
  end_date <- planting_date + days(max_growing_days - 1)
  
  sim_period <- weather %>%
    filter(date >= planting_date, date <= end_date) %>%
    arrange(date)
  
  if (nrow(sim_period) == 0) {
    return(data.frame(
      crop_name = crop_name, crop_type = crop_type, row_role = "harvest_year", year = yr,
      planting_date = planting_date, maturity_date = as.Date(NA), harvest_date = as.Date(NA),
      season_length_days = NA_integer_, accumulated_tt = 0, required_tt = required_tt,
      maturity_fraction = 0, status = "failed_after_planting", forced_harvest = FALSE,
      regrowth_started = FALSE,
      vernalization_satisfied = ifelse(crop_type == "winter" && vernalization_required, FALSE, NA),
      vernalization_days = NA_real_
    ))
  }
  
  sim_period$daily_tt <- 0
  sim_period$cum_tt <- 0
  vernal_days <- 0
  vern_sat <- ifelse(crop_type == "winter" && vernalization_required, FALSE, TRUE)
  regrowth_started <- FALSE
  cum_tt <- 0
  
  for (i in seq_len(nrow(sim_period))) {
    tmean_i <- sim_period$tmean[i]
    
    if (crop_type == "winter" && vernalization_required && !vern_sat) {
      if (tmean_i >= vernalization_temp_min && tmean_i <= vernalization_temp_max) vernal_days <- vernal_days + 1
      if (vernal_days >= vernalization_days_required) vern_sat <- TRUE
    }
    
    if (crop_type == "winter" && vern_sat && !regrowth_started && tmean_i >= spring_regrowth_temp) {
      regrowth_started <- TRUE
    }
    
    if (crop_type == "winter" && !vern_sat) {
      tt_i <- 0
    } else {
      tt_i <- calc_daily_tt(sim_period$tmin[i], sim_period$tmax[i], t_base, t_opt, t_max_cut, mode = tt_mode)
      if (crop_type == "winter" && tmean_i <= winter_dormancy_temp) tt_i <- 0
    }
    
    cum_tt <- cum_tt + tt_i
    sim_period$daily_tt[i] <- tt_i
    sim_period$cum_tt[i] <- cum_tt
  }
  
  accumulated_tt <- dplyr::last(sim_period$cum_tt)
  maturity_fraction <- accumulated_tt / required_tt
  mature_idx <- which(sim_period$cum_tt >= required_tt)
  
  if (length(mature_idx) > 0) {
    maturity_date <- sim_period$date[min(mature_idx)]
    harvest_date <- maturity_date
    accumulated_tt <- sim_period$cum_tt[min(mature_idx)]
    maturity_fraction <- accumulated_tt / required_tt
    
    return(data.frame(
      crop_name = crop_name, crop_type = crop_type, row_role = "harvest_year", year = yr,
      planting_date = planting_date, maturity_date = maturity_date, harvest_date = harvest_date,
      season_length_days = as.integer(harvest_date - planting_date) + 1,
      accumulated_tt = accumulated_tt, required_tt = required_tt, maturity_fraction = maturity_fraction,
      status = "mature", forced_harvest = FALSE, regrowth_started = regrowth_started,
      vernalization_satisfied = ifelse(crop_type == "winter" && vernalization_required, vern_sat, NA),
      vernalization_days = ifelse(crop_type == "winter" && vernalization_required, vernal_days, NA)
    ))
  }
  
  if (forced_harvest_allowed && maturity_fraction >= min_fraction_tt_for_forced_harvest) {
    harvest_date <- dplyr::last(sim_period$date)
    return(data.frame(
      crop_name = crop_name, crop_type = crop_type, row_role = "harvest_year", year = yr,
      planting_date = planting_date, maturity_date = as.Date(NA), harvest_date = harvest_date,
      season_length_days = as.integer(harvest_date - planting_date) + 1,
      accumulated_tt = accumulated_tt, required_tt = required_tt, maturity_fraction = maturity_fraction,
      status = "forced_harvest_immature", forced_harvest = TRUE, regrowth_started = regrowth_started,
      vernalization_satisfied = ifelse(crop_type == "winter" && vernalization_required, vern_sat, NA),
      vernalization_days = ifelse(crop_type == "winter" && vernalization_required, vernal_days, NA)
    ))
  }
  
  data.frame(
    crop_name = crop_name, crop_type = crop_type, row_role = "harvest_year", year = yr,
    planting_date = planting_date, maturity_date = as.Date(NA), harvest_date = as.Date(NA),
    season_length_days = as.integer(dplyr::last(sim_period$date) - planting_date) + 1,
    accumulated_tt = accumulated_tt, required_tt = required_tt, maturity_fraction = maturity_fraction,
    status = ifelse(crop_type == "winter" && vernalization_required && !vern_sat, "insufficient_vernalization", "failed_to_mature"),
    forced_harvest = FALSE, regrowth_started = regrowth_started,
    vernalization_satisfied = ifelse(crop_type == "winter" && vernalization_required, vern_sat, NA),
    vernalization_days = ifelse(crop_type == "winter" && vernalization_required, vernal_days, NA)
  )
}

run_simulation <- function(weather,
                           crop_name,
                           required_tt,
                           max_growing_days,
                           t_base,
                           t_opt = NA,
                           t_max_cut = NA,
                           tt_mode = "simple",
                           crop_type = c("summer", "winter"),
                           min_mean_temp_plant = 0,
                           forced_harvest_allowed = TRUE,
                           min_fraction_tt_for_forced_harvest = 0,
                           winter_dormancy_temp = 0,
                           vernalization_required = FALSE,
                           vernalization_temp_min = 0,
                           vernalization_temp_max = 10,
                           vernalization_days_required = 30,
                           spring_regrowth_temp = 5,
                           winter_plant_temp_min = 5,
                           winter_plant_temp_max = 15) {
  crop_type <- match.arg(crop_type)
  yrs <- sort(unique(weather$year))
  
  if (crop_type == "summer") {
    out <- bind_rows(lapply(yrs, function(yr) {
      simulate_one_year(
        weather, yr, required_tt, max_growing_days,
        t_base, t_opt, t_max_cut, tt_mode, crop_type,
        min_mean_temp_plant, forced_harvest_allowed, min_fraction_tt_for_forced_harvest,
        winter_dormancy_temp, vernalization_required, vernalization_temp_min,
        vernalization_temp_max, vernalization_days_required, spring_regrowth_temp,
        winter_plant_temp_min, winter_plant_temp_max, crop_name
      )
    }))
  } else {
    first_year <- min(yrs)
    last_year <- max(yrs)
    
    planting_rows <- bind_rows(lapply(yrs[yrs < last_year], function(yr) {
      windows <- get_internal_windows(crop_type)
      weather_window <- weather %>%
        filter(date >= clamp_year_day(yr, windows$earliest_planting_mmdd),
               date <= clamp_year_day(yr, windows$latest_planting_mmdd))
      
      plant_result <- find_planting_date(
        weather_window,
        clamp_year_day(yr, windows$earliest_planting_mmdd),
        clamp_year_day(yr, windows$latest_planting_mmdd),
        crop_type,
        min_mean_temp_plant,
        winter_plant_temp_min,
        winter_plant_temp_max
      )
      
      data.frame(
        crop_name = crop_name, crop_type = crop_type, row_role = "planting_year", year = yr,
        planting_date = if (plant_result$found) plant_result$planting_date else as.Date(NA),
        maturity_date = as.Date(NA), harvest_date = as.Date(NA), season_length_days = NA_integer_,
        accumulated_tt = NA_real_, required_tt = required_tt, maturity_fraction = NA_real_,
        status = if (plant_result$found) "planted_for_next_year" else "not_planted",
        forced_harvest = FALSE, regrowth_started = FALSE,
        vernalization_satisfied = NA, vernalization_days = NA_real_
      )
    }))
    
    harvest_rows <- bind_rows(lapply(yrs[yrs > first_year], function(yr) {
      full_res <- simulate_one_year(
        weather, yr, required_tt, max_growing_days,
        t_base, t_opt, t_max_cut, tt_mode, crop_type,
        min_mean_temp_plant, forced_harvest_allowed, min_fraction_tt_for_forced_harvest,
        winter_dormancy_temp, vernalization_required, vernalization_temp_min,
        vernalization_temp_max, vernalization_days_required, spring_regrowth_temp,
        winter_plant_temp_min, winter_plant_temp_max, crop_name
      )
      full_res$planting_date <- as.Date(NA)
      full_res$row_role <- "harvest_year"
      full_res
    }))
    
    out <- bind_rows(planting_rows, harvest_rows) %>% arrange(year, row_role)
  }
  
  out
}

build_events_plot_data <- function(sim_results) {
  if (nrow(sim_results) == 0) return(list(events_plot = NULL, bx = NULL, crop_lab = NULL))
  
  x <- sim_results
  if (!"dataset" %in% names(x)) x$dataset <- "Historical"
  if (!"scenario" %in% names(x)) x$scenario <- NA_character_
  if (!"model" %in% names(x)) x$model <- NA_character_
  
  x <- x %>%
    mutate(
      Scenario = case_when(
        dataset == "Historical" ~ "Historical",
        !is.na(scenario) & !is.na(model) ~ paste0(scenario, " | ", model),
        !is.na(scenario) ~ scenario,
        !is.na(model) ~ model,
        TRUE ~ dataset
      )
    )
  
  events <- bind_rows(
    x %>% filter(!is.na(planting_date)) %>% transmute(crop = crop_name, Scenario, operation = factor("PLANT", levels = c("PLANT", "HARV/KILL")), date = planting_date, doy = yday(planting_date)),
    x %>% filter(!is.na(harvest_date)) %>% transmute(crop = crop_name, Scenario, operation = factor("HARV/KILL", levels = c("PLANT", "HARV/KILL")), date = harvest_date, doy = yday(harvest_date))
  )
  
  if (nrow(events) == 0) return(list(events_plot = NULL, bx = NULL, crop_lab = NULL))
  
  events <- events %>%
    mutate(
      crop = factor(crop, levels = unique(crop)),
      Scenario = factor(Scenario, levels = unique(Scenario)),
      operation = factor(operation, levels = c("PLANT", "HARV/KILL"))
    )
  
  events_clean <- events %>%
    group_by(crop, Scenario, operation) %>%
    mutate(
      q1 = quantile(doy, 0.25, na.rm = TRUE),
      q3 = quantile(doy, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower = q1 - 1.5 * iqr,
      upper = q3 + 1.5 * iqr
    ) %>%
    ungroup() %>%
    filter(doy >= lower, doy <= upper)
  
  events_plot <- events_clean %>%
    mutate(
      doy_plot = doy,
      Scenario = droplevels(Scenario),
      Scenario_y = as.numeric(Scenario),
      op_off = if_else(operation == "PLANT", 0.18, -0.18),
      y_ridge = Scenario_y + op_off
    )
  
  bx <- events_plot %>%
    group_by(crop, Scenario, Scenario_y, operation, op_off) %>%
    summarise(
      q1 = quantile(doy_plot, 0.25, na.rm = TRUE),
      q2 = quantile(doy_plot, 0.50, na.rm = TRUE),
      q3 = quantile(doy_plot, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      low = min(doy_plot[doy_plot >= (q1 - 1.5 * iqr)], na.rm = TRUE),
      up = max(doy_plot[doy_plot <= (q3 + 1.5 * iqr)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_box = Scenario_y + op_off - 0.12)
  
  crop_lab <- setNames(as.character(unique(events_plot$crop)), unique(events_plot$crop))
  
  list(events_plot = events_plot, bx = bx, crop_lab = crop_lab)
}

plot_growing_season_ridges <- function(sim_results) {
  obj <- build_events_plot_data(sim_results)
  events_plot <- obj$events_plot
  bx <- obj$bx
  crop_lab <- obj$crop_lab
  
  if (is.null(events_plot) || nrow(events_plot) == 0) {
    return(ggplot() + theme_void() + annotate("text", x = 0, y = 0, label = "No valid event dates to plot."))
  }
  
  event_cols <- brewer.pal(3, "Pastel1")[1:2]
  names(event_cols) <- c("PLANT", "HARV/KILL")
  
  theme_pub <- theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey85"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.box = "horizontal",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
  
  month_dates <- ymd(paste0("2001-", sprintf("%02d", 1:12), "-01"))
  month_breaks <- yday(month_dates)
  month_labels <- month(month_dates, label = TRUE, abbr = TRUE)
  
  box_h <- 0.055
  cap_h <- 0.030
  
  ggplot(
    events_plot,
    aes(x = doy_plot, y = y_ridge, fill = operation,
        group = interaction(Scenario, operation, crop))
  ) +
    ggridges::geom_density_ridges(
      alpha = 0.88, scale = 1.12, rel_min_height = 0.01,
      linewidth = 0.35, colour = "grey35", bandwidth = 4
    ) +
    geom_segment(data = bx, aes(x = low, xend = q1, y = y_box, yend = y_box),
                 inherit.aes = FALSE, linewidth = 0.35, colour = "grey35") +
    geom_segment(data = bx, aes(x = q3, xend = up, y = y_box, yend = y_box),
                 inherit.aes = FALSE, linewidth = 0.35, colour = "grey35") +
    geom_segment(data = bx, aes(x = low, xend = low, y = y_box - cap_h, yend = y_box + cap_h),
                 inherit.aes = FALSE, linewidth = 0.35, colour = "grey35") +
    geom_segment(data = bx, aes(x = up, xend = up, y = y_box - cap_h, yend = y_box + cap_h),
                 inherit.aes = FALSE, linewidth = 0.35, colour = "grey35") +
    geom_rect(
      data = bx,
      aes(xmin = q1, xmax = q3, ymin = y_box - box_h, ymax = y_box + box_h, fill = operation),
      inherit.aes = FALSE, alpha = 0.95, colour = "grey40", linewidth = 0.35
    ) +
    geom_segment(data = bx, aes(x = q2, xend = q2, y = y_box - box_h, yend = y_box + box_h),
                 inherit.aes = FALSE, linewidth = 0.45, colour = "grey20") +
    facet_wrap(~ crop, ncol = 1, labeller = labeller(crop = crop_lab)) +
    scale_x_continuous(
      "Calendar month",
      breaks = month_breaks, labels = month_labels, limits = c(1, 365),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      "Dataset / scenario",
      breaks = sort(unique(events_plot$Scenario_y)),
      labels = levels(events_plot$Scenario),
      expand = expansion(mult = c(0.10, 0.12))
    ) +
    scale_fill_manual(
      values = event_cols, breaks = c("PLANT", "HARV/KILL"),
      labels = c("Planting", "Harvest"), name = "Field operation"
    ) +
    labs(
      title = "Growing season timing across historical and climate scenarios",
      subtitle = "Density ridges show the distribution of planting and harvest dates; inline boxplots summarize median, interquartile range, and whiskers for each operation."
    ) +
    coord_cartesian(clip = "off") +
    theme_pub
}

compare_summary <- function(combined_results) {
  x <- combined_results
  if (nrow(x) == 0) return(data.frame())
  
  if (!"dataset" %in% names(x)) x$dataset <- "Historical"
  if (!"scenario" %in% names(x)) x$scenario <- NA_character_
  if (!"model" %in% names(x)) x$model <- NA_character_
  
  x %>%
    mutate(
      group_label = case_when(
        dataset == "Historical" ~ "Historical",
        !is.na(scenario) & !is.na(model) ~ paste0(scenario, " | ", model),
        !is.na(scenario) ~ scenario,
        !is.na(model) ~ model,
        TRUE ~ dataset
      )
    ) %>%
    filter(row_role == "harvest_year") %>%
    group_by(group_label) %>%
    summarise(
      n_years = n(),
      mature_pct = round(mean(status == "mature", na.rm = TRUE) * 100, 1),
      forced_harvest_pct = round(mean(status == "forced_harvest_immature", na.rm = TRUE) * 100, 1),
      failed_pct = round(mean(status %in% c("failed_to_mature", "insufficient_vernalization"), na.rm = TRUE) * 100, 1),
      median_harvest_doy = round(median(yday(harvest_date), na.rm = TRUE), 1),
      median_maturity_doy = round(median(yday(maturity_date), na.rm = TRUE), 1),
      mean_maturity_fraction = round(mean(maturity_fraction, na.rm = TRUE), 3),
      mean_planting_to_harvest_days = round(mean(season_length_days, na.rm = TRUE), 1),
      median_planting_to_harvest_days = round(median(season_length_days, na.rm = TRUE), 1),
      .groups = "drop"
    )
}

make_summary_text <- function(combined_results) {
  x <- compare_summary(combined_results)
  if (nrow(x) == 0) return("No results available.")
  
  hist <- x %>% filter(group_label == "Historical")
  fut <- x %>% filter(group_label != "Historical")
  
  if (nrow(hist) == 0) return("No historical summary available.")
  if (nrow(fut) == 0) {
    return(paste0(
      "Historical simulation only. Mature years: ", hist$mature_pct[1], "%; ",
      "median harvest day-of-year: ", hist$median_harvest_doy[1], "; ",
      "median season length: ", hist$median_planting_to_harvest_days[1], " days."
    ))
  }
  
  best_early <- fut %>% arrange(median_harvest_doy) %>% slice(1)
  best_short <- fut %>% arrange(median_planting_to_harvest_days) %>% slice(1)
  
  paste0(
    "Historical baseline: ", hist$mature_pct[1], "% mature years, median harvest DOY ",
    hist$median_harvest_doy[1], ", median season length ",
    hist$median_planting_to_harvest_days[1], " days. ",
    "Across climate scenarios, the earliest median harvest occurs under ",
    best_early$group_label[1], " (DOY ", best_early$median_harvest_doy[1], "). ",
    "The shortest median planting-to-harvest duration occurs under ",
    best_short$group_label[1], " (", best_short$median_planting_to_harvest_days[1], " days)."
  )
}


prepare_temperature_plot_data <- function(hist_weather, climate_weather = NULL) {
  hw <- hist_weather %>%
    mutate(
      dataset = "Historical",
      group_label = "Historical",
      month = month(date, label = TRUE, abbr = TRUE),
      month_num = month(date)
    )
  
  if (is.null(climate_weather)) {
    cw <- NULL
  } else {
    cw <- climate_weather
    if (!"scenario" %in% names(cw)) cw$scenario <- NA_character_
    if (!"model" %in% names(cw)) cw$model <- NA_character_
    cw <- cw %>%
      mutate(
        dataset = "Climate",
        group_label = case_when(
          !is.na(scenario) & !is.na(model) ~ paste0(scenario, " | ", model),
          !is.na(scenario) ~ scenario,
          !is.na(model) ~ model,
          TRUE ~ "Climate"
        ),
        month = month(date, label = TRUE, abbr = TRUE),
        month_num = month(date)
      )
  }
  
  bind_rows(hw, cw)
}

plot_monthly_temperature_cycle <- function(hist_weather, climate_weather = NULL) {
  df <- prepare_temperature_plot_data(hist_weather, climate_weather)
  
  monthly <- df %>%
    group_by(group_label, dataset, month, month_num) %>%
    summarise(
      mean_tmean = mean(tmean, na.rm = TRUE),
      q10 = quantile(tmean, 0.10, na.rm = TRUE),
      q90 = quantile(tmean, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(monthly, aes(x = month_num, y = mean_tmean, group = group_label, color = group_label, fill = group_label)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.12, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 1.8) +
    scale_x_continuous(
      breaks = 1:12,
      labels = month.abb,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Seasonal temperature cycle",
      subtitle = "Monthly mean temperature with 10th–90th percentile envelope",
      x = "Month",
      y = "Mean daily temperature (°C)",
      color = "Scenario",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

plot_annual_temperature_boxplot <- function(hist_weather, climate_weather = NULL) {
  df <- prepare_temperature_plot_data(hist_weather, climate_weather)
  
  annual <- df %>%
    group_by(group_label, year) %>%
    summarise(
      annual_mean_tmean = mean(tmean, na.rm = TRUE),
      .groups = "drop"
    )
  ggplot(annual, aes(
    x = reorder(group_label, annual_mean_tmean, FUN = median),
    y = annual_mean_tmean,
    fill = group_label
  )) +
    geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
    coord_flip() +
    scale_fill_hue(name = "Scenario") +
    guides(fill = guide_legend(ncol = 1)) +
    labs(
      x = NULL,
      y = "Annual mean temperature (°C)",
      title = "Annual temperature distribution",
      subtitle = "Each box summarizes annual mean temperature across years"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )

}

plot_growing_season_heatmap <- function(hist_weather, climate_weather = NULL) {
  df <- prepare_temperature_plot_data(hist_weather, climate_weather)
  
  monthly <- df %>%
    group_by(group_label, month, month_num) %>%
    summarise(mean_tmean = mean(tmean, na.rm = TRUE), .groups = "drop")
  
  ggplot(monthly, aes(x = month_num, y = reorder(group_label, mean_tmean, FUN = mean), fill = mean_tmean)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +
    labs(
      x = NULL,
      y = NULL,
      fill = "°C",
      title = "Monthly temperature heatmap",
      subtitle = "Fast comparison of historical and projected temperature regime"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}



compute_tt_diagnostics <- function(weather, group_label,
                                   t_base, t_opt = NA, t_max_cut = NA,
                                   tt_mode = "simple") {
  weather %>%
    mutate(
      group_label = group_label,
      daily_tt = calc_daily_tt(tmin, tmax, t_base, t_opt, t_max_cut, mode = tt_mode),
      annual_tt = ave(daily_tt, year, FUN = sum),
      month = month(date, label = TRUE, abbr = TRUE),
      month_num = month(date)
    )
}

build_tt_diagnostics_data <- function(hist_weather, climate_weather = NULL,
                                      t_base, t_opt = NA, t_max_cut = NA,
                                      tt_mode = "simple") {
  hist_df <- compute_tt_diagnostics(hist_weather, "Historical", t_base, t_opt, t_max_cut, tt_mode)
  
  if (is.null(climate_weather)) return(hist_df)
  
  cw <- climate_weather
  if (!"scenario" %in% names(cw)) cw$scenario <- NA_character_
  if (!"model" %in% names(cw)) cw$model <- NA_character_
  
  pieces <- cw %>%
    mutate(
      group_label = case_when(
        !is.na(scenario) & !is.na(model) ~ paste0(scenario, " | ", model),
        !is.na(scenario) ~ scenario,
        !is.na(model) ~ model,
        TRUE ~ "Climate"
      )
    ) %>%
    group_by(group_label) %>%
    group_split()
  
  fut_df <- dplyr::bind_rows(lapply(pieces, function(d) {
    compute_tt_diagnostics(d, unique(d$group_label), t_base, t_opt, t_max_cut, tt_mode)
  }))
  
  bind_rows(hist_df, fut_df)
}

plot_monthly_tt_cycle <- function(hist_weather, climate_weather = NULL,
                                  t_base, t_opt = NA, t_max_cut = NA,
                                  tt_mode = "simple") {
  df <- build_tt_diagnostics_data(hist_weather, climate_weather, t_base, t_opt, t_max_cut, tt_mode)
  
  monthly <- df %>%
    group_by(group_label, month, month_num) %>%
    summarise(
      mean_daily_tt = mean(daily_tt, na.rm = TRUE),
      q10 = quantile(daily_tt, 0.10, na.rm = TRUE),
      q90 = quantile(daily_tt, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(monthly, aes(x = month_num, y = mean_daily_tt, group = group_label, color = group_label, fill = group_label)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.12, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 1.7) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(mult = c(0.01, 0.02))) +
    labs(
      title = "Seasonal thermal-time cycle",
      subtitle = "Monthly mean daily thermal time with 10th–90th percentile envelope",
      x = "Month",
      y = "Daily thermal time (TT units day⁻¹)",
      color = "Scenario",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

plot_annual_tt_boxplot <- function(hist_weather, climate_weather = NULL,
                                   t_base, t_opt = NA, t_max_cut = NA,
                                   tt_mode = "simple") {
  df <- build_tt_diagnostics_data(hist_weather, climate_weather, t_base, t_opt, t_max_cut, tt_mode)
  
  annual <- df %>%
    group_by(group_label, year) %>%
    summarise(annual_tt = sum(daily_tt, na.rm = TRUE), .groups = "drop")
  
  ggplot(annual, aes(x = reorder(group_label, annual_tt, FUN = median), y = annual_tt, fill = group_label)) +
    geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
    coord_flip() +
    labs(
      title = "Annual accumulated thermal time",
      subtitle = "Distribution of yearly accumulated thermal time across scenarios",
      x = NULL,
      y = "Annual TT sum",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

plot_tmin_tmax_ribbon <- function(hist_weather, climate_weather = NULL) {
  df <- prepare_temperature_plot_data(hist_weather, climate_weather)
  
  monthly <- df %>%
    group_by(group_label, month, month_num) %>%
    summarise(
      mean_tmin = mean(tmin, na.rm = TRUE),
      mean_tmax = mean(tmax, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(monthly, aes(x = month_num, group = group_label, color = group_label, fill = group_label)) +
    geom_ribbon(aes(ymin = mean_tmin, ymax = mean_tmax), alpha = 0.12, linewidth = 0, show.legend = FALSE) +
    geom_line(aes(y = mean_tmin), linewidth = 0.75, linetype = "dashed") +
    geom_line(aes(y = mean_tmax), linewidth = 0.95) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(mult = c(0.01, 0.02))) +
    labs(
      title = "Monthly Tmin–Tmax envelope",
      subtitle = "Dashed line = mean Tmin, solid line = mean Tmax",
      x = "Month",
      y = "Temperature (°C)",
      color = "Scenario",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

ui <- fluidPage(
  titlePanel("PhenoFlex - Historical + Climate Change app"),
  sidebarLayout(
    sidebarPanel(
      fileInput("weather_file", "Upload historical weather CSV", accept = ".csv"),
      fileInput("climate_file", "Upload climate-change CSV (optional)", accept = ".csv"),
      tags$small("Required columns: date, tmin, tmax. Optional climate metadata: scenario, model, period, station."),
      hr(),
      
      h4("Detected baseline period"),
      verbatimTextOutput("baseline_years_text"),
      hr(),
      
      h4("Crop settings"),
      textInput("crop_name", "Crop name", value = "Maize"),
      selectInput("crop_type", "Crop type", choices = c("summer", "winter"), selected = "summer"),
      
      h4("Crop thermal settings"),
      numericInput("days_to_maturity", "Days to maturity", value = 140, min = 1),
      numericInput("t_base", "Base temperature (°C)", value = 8),
      numericInput("t_opt", "Optimum temperature (°C)", value = 25),
      numericInput("t_max_cut", "Upper temperature cutoff (°C)", value = 35),
      selectInput("tt_mode", "Thermal time method", choices = c("simple", "capped", "triangular"), selected = "triangular"),
      textInput("baseline_planting_mmdd", "Reference planting date (MM-DD)", value = "04-15"),
      numericInput("max_growing_days", "Maximum growing days", value = 220, min = 1),
      hr(),
      
      h4("Planting thresholds"),
      numericInput("min_mean_temp_plant", "Summer minimum mean temperature for planting (°C)", value = 8),
      numericInput("winter_plant_temp_min", "Winter crop planting mean temp min (°C)", value = 5),
      numericInput("winter_plant_temp_max", "Winter crop planting mean temp max (°C)", value = 15),
      hr(),
      
      conditionalPanel(
        condition = "input.crop_type == 'winter'",
        h4("Winter crop settings"),
        numericInput("winter_dormancy_temp", "Winter dormancy temperature (°C)", value = 0),
        checkboxInput("vernalization_required", "Require vernalization", value = TRUE),
        numericInput("vernalization_temp_min", "Vernalization min temperature (°C)", value = 0),
        numericInput("vernalization_temp_max", "Vernalization max temperature (°C)", value = 10),
        numericInput("vernalization_days_required", "Vernalization days required", value = 45),
        numericInput("spring_regrowth_temp", "Spring regrowth temperature (°C)", value = 5)
      ),
      hr(),
      
      h4("Harvest rule"),
      checkboxInput("forced_harvest_allowed", "Allow forced harvest at max growing days", value = TRUE),
      numericInput("min_fraction_tt_for_forced_harvest", "Minimum maturity fraction for forced harvest", value = 0.8, min = 0, max = 1, step = 0.05),
      hr(),
      
      actionButton("run_model", "Run simulation")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Overview",
                 h4("Run settings"),
                 p("Run settings describe the crop, thermal-development assumptions, and baseline reference used by the model before running the simulations. These settings control how daily temperature is converted into crop development, how the required thermal time is estimated from the historical baseline, and how the crop is then simulated under historical and climate-change conditions."),
                 tags$ul(
                   tags$li(tags$b("Crop name:"), " The crop being simulated."),
                   tags$li(tags$b("Crop type:"), " Indicates whether the crop is treated as a summer crop or a winter crop. Winter crops include autumn planting, overwintering, vernalization, and spring regrowth logic."),
                   tags$li(tags$b("Days to maturity:"), " The reference number of calendar days used to estimate the crop’s required thermal time from the historical baseline dataset."),
                   tags$li(tags$b("Maximum growing days:"), " The maximum allowed crop season length in the simulation. If the crop does not reach maturity before this limit, the model may classify it as failed or force harvest if that option is enabled."),
                   tags$li(tags$b("Base temperature (°C):"), " The lower threshold for thermal development. Temperatures below this value do not contribute to crop development."),
                   tags$li(tags$b("Optimum temperature (°C):"), " The temperature at which thermal development is assumed to be most effective."),
                   tags$li(tags$b("Upper temperature cutoff (°C):"), " The temperature above which thermal development no longer increases effectively in the selected thermal-time response method."),
                   tags$li(tags$b("Thermal time method:"), " The mathematical rule used to convert daily temperature into daily crop development."),
                   tags$li(tags$b("Baseline years:"), " The historical period used to estimate the crop’s required thermal time."),
                   tags$li(tags$b("Reference planting date:"), " The planting date used in the baseline historical period to estimate required thermal time from the selected number of maturity days.")
                 ),
                 tableOutput("run_parameters_table"),
                 h4("Estimated required thermal time"),
                 verbatimTextOutput("required_tt_text"),
                 h4("Comparison summary"),
                 p("Comparison summary shows how the crop behaves under the historical baseline and under each climate scenario after applying the same crop parameters and the same required thermal time estimated from the historical period. The table summarizes how often the crop matures successfully, how early or late harvest occurs, and how the total growing-season length changes across scenarios."),
                 tags$ul(
                   tags$li(tags$b("group_label:"), " The dataset or climate scenario being summarized."),
                   tags$li(tags$b("n_years:"), " The number of simulated years included in that summary group."),
                   tags$li(tags$b("mature_pct:"), " Percentage of years in which the crop reached full maturity before the end of the allowed growing period."),
                   tags$li(tags$b("forced_harvest_pct:"), " Percentage of years in which the crop did not fully mature but was still harvested because it had reached at least the minimum maturity fraction required for forced harvest."),
                   tags$li(tags$b("failed_pct:"), " Percentage of years in which the crop did not mature and was not eligible for forced harvest."),
                   tags$li(tags$b("median_harvest_doy:"), " Median day of year for harvest across all simulated years in that group. Lower values mean earlier harvest."),
                   tags$li(tags$b("median_maturity_doy:"), " Median day of year on which the crop reached maturity. Lower values mean earlier maturity."),
                   tags$li(tags$b("mean_maturity_fraction:"), " Average fraction of the required thermal time reached by the crop. Values around 1.00 mean the crop usually reached full maturity."),
                   tags$li(tags$b("mean_planting_to_harvest_days:"), " Average number of days from planting to harvest across simulated years."),
                   tags$li(tags$b("median_planting_to_harvest_days:"), " Median number of days from planting to harvest across simulated years.")
                 ),
                 tableOutput("comparison_summary_table")),
        tabPanel("Combined results", DTOutput("results_table")),
        tabPanel("Plots",
                 h4("Interpretation"),
                 verbatimTextOutput("summary_text"),
                 plotOutput("plot_ridge_boxinline", height = 700),
                 hr(),
                 h4("Climate diagnostics"),
                 plotOutput("plot_monthly_temperature_cycle", height = 500),
                 plotOutput("plot_tmin_tmax_ribbon", height = 500),
                 plotOutput("plot_annual_temperature_boxplot", height = 430),
                 plotOutput("plot_monthly_tt_cycle", height = 500),
                 plotOutput("plot_annual_tt_boxplot", height = 430),
                 plotOutput("plot_temperature_heatmap", height = 420)),
        tabPanel("Historical weather",
                 DTOutput("weather_preview"),
                 verbatimTextOutput("weather_info")),
        tabPanel("Projected weather",
                 DTOutput("climate_preview"),
                 verbatimTextOutput("climate_info"))
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$crop_type, {
    p <- default_parameters(input$crop_type)
    updateTextInput(session, "crop_name", value = p$crop_name)
    updateNumericInput(session, "days_to_maturity", value = p$days_to_maturity)
    updateNumericInput(session, "t_base", value = p$t_base)
    updateNumericInput(session, "t_opt", value = p$t_opt)
    updateNumericInput(session, "t_max_cut", value = p$t_max_cut)
    updateTextInput(session, "baseline_planting_mmdd", value = p$baseline_planting_mmdd)
    updateNumericInput(session, "max_growing_days", value = p$max_growing_days)
    updateNumericInput(session, "min_mean_temp_plant", value = p$min_mean_temp_plant)
    updateNumericInput(session, "winter_plant_temp_min", value = p$winter_plant_temp_min)
    updateNumericInput(session, "winter_plant_temp_max", value = p$winter_plant_temp_max)
  }, ignoreInit = TRUE)
  
  weather_data <- reactive({
    req(input$weather_file)
    prepare_weather(read_csv(input$weather_file$datapath, show_col_types = FALSE))
  })
  
  climate_data <- reactive({
    req(input$climate_file)
    prepare_weather(read_csv(input$climate_file$datapath, show_col_types = FALSE))
  })
  
  baseline_years_auto <- reactive({
    req(weather_data())
    seq(min(weather_data()$year), max(weather_data()$year))
  })
  
  output$baseline_years_text <- renderPrint({
    req(weather_data())
    yrs <- baseline_years_auto()
    cat("Baseline years:", min(yrs), "-", max(yrs), "\n")
  })
  
  required_tt_result <- eventReactive(input$run_model, {
    estimate_required_tt(
      weather = weather_data(),
      baseline_years = baseline_years_auto(),
      planting_mmdd = input$baseline_planting_mmdd,
      days_to_maturity = input$days_to_maturity,
      t_base = input$t_base,
      t_opt = input$t_opt,
      t_max_cut = input$t_max_cut,
      tt_mode = input$tt_mode,
      crop_type = input$crop_type,
      winter_dormancy_temp = ifelse(input$crop_type == "winter", input$winter_dormancy_temp, 0),
      vernalization_required = ifelse(input$crop_type == "winter", input$vernalization_required, FALSE),
      vernalization_temp_min = ifelse(input$crop_type == "winter", input$vernalization_temp_min, 0),
      vernalization_temp_max = ifelse(input$crop_type == "winter", input$vernalization_temp_max, 10),
      vernalization_days_required = ifelse(input$crop_type == "winter", input$vernalization_days_required, 0),
      spring_regrowth_temp = ifelse(input$crop_type == "winter", input$spring_regrowth_temp, 5)
    )
  })
  
  historical_results <- eventReactive(input$run_model, {
    run_simulation(
      weather = weather_data(),
      crop_name = input$crop_name,
      required_tt = required_tt_result()$required_tt,
      max_growing_days = input$max_growing_days,
      t_base = input$t_base,
      t_opt = input$t_opt,
      t_max_cut = input$t_max_cut,
      tt_mode = input$tt_mode,
      crop_type = input$crop_type,
      min_mean_temp_plant = input$min_mean_temp_plant,
      forced_harvest_allowed = input$forced_harvest_allowed,
      min_fraction_tt_for_forced_harvest = input$min_fraction_tt_for_forced_harvest,
      winter_dormancy_temp = ifelse(input$crop_type == "winter", input$winter_dormancy_temp, 0),
      vernalization_required = ifelse(input$crop_type == "winter", input$vernalization_required, FALSE),
      vernalization_temp_min = ifelse(input$crop_type == "winter", input$vernalization_temp_min, 0),
      vernalization_temp_max = ifelse(input$crop_type == "winter", input$vernalization_temp_max, 10),
      vernalization_days_required = ifelse(input$crop_type == "winter", input$vernalization_days_required, 0),
      spring_regrowth_temp = ifelse(input$crop_type == "winter", input$spring_regrowth_temp, 5),
      winter_plant_temp_min = input$winter_plant_temp_min,
      winter_plant_temp_max = input$winter_plant_temp_max
    ) %>%
      mutate(dataset = "Historical", scenario = NA_character_, model = NA_character_)
  })
  
  future_results <- eventReactive(input$run_model, {
    req(input$climate_file)
    cd <- climate_data()
    
    grouping_cols <- intersect(c("scenario", "model", "period", "station"), names(cd))
    
    if (length(grouping_cols) == 0) {
      run_simulation(
        weather = cd,
        crop_name = input$crop_name,
        required_tt = required_tt_result()$required_tt,
        max_growing_days = input$max_growing_days,
        t_base = input$t_base,
        t_opt = input$t_opt,
        t_max_cut = input$t_max_cut,
        tt_mode = input$tt_mode,
        crop_type = input$crop_type,
        min_mean_temp_plant = input$min_mean_temp_plant,
        forced_harvest_allowed = input$forced_harvest_allowed,
        min_fraction_tt_for_forced_harvest = input$min_fraction_tt_for_forced_harvest,
        winter_dormancy_temp = ifelse(input$crop_type == "winter", input$winter_dormancy_temp, 0),
        vernalization_required = ifelse(input$crop_type == "winter", input$vernalization_required, FALSE),
        vernalization_temp_min = ifelse(input$crop_type == "winter", input$vernalization_temp_min, 0),
        vernalization_temp_max = ifelse(input$crop_type == "winter", input$vernalization_temp_max, 10),
        vernalization_days_required = ifelse(input$crop_type == "winter", input$vernalization_days_required, 0),
        spring_regrowth_temp = ifelse(input$crop_type == "winter", input$spring_regrowth_temp, 5),
        winter_plant_temp_min = input$winter_plant_temp_min,
        winter_plant_temp_max = input$winter_plant_temp_max
      ) %>%
        mutate(dataset = "Climate", scenario = NA_character_, model = NA_character_)
    } else {
      cd %>%
        group_by(across(all_of(grouping_cols))) %>%
        group_split() %>%
        lapply(function(df_group) {
          meta <- df_group[1, grouping_cols, drop = FALSE]
          
          res <- run_simulation(
            weather = df_group,
            crop_name = input$crop_name,
            required_tt = required_tt_result()$required_tt,
            max_growing_days = input$max_growing_days,
            t_base = input$t_base,
            t_opt = input$t_opt,
            t_max_cut = input$t_max_cut,
            tt_mode = input$tt_mode,
            crop_type = input$crop_type,
            min_mean_temp_plant = input$min_mean_temp_plant,
            forced_harvest_allowed = input$forced_harvest_allowed,
            min_fraction_tt_for_forced_harvest = input$min_fraction_tt_for_forced_harvest,
            winter_dormancy_temp = ifelse(input$crop_type == "winter", input$winter_dormancy_temp, 0),
            vernalization_required = ifelse(input$crop_type == "winter", input$vernalization_required, FALSE),
            vernalization_temp_min = ifelse(input$crop_type == "winter", input$vernalization_temp_min, 0),
            vernalization_temp_max = ifelse(input$crop_type == "winter", input$vernalization_temp_max, 10),
            vernalization_days_required = ifelse(input$crop_type == "winter", input$vernalization_days_required, 0),
            spring_regrowth_temp = ifelse(input$crop_type == "winter", input$spring_regrowth_temp, 5),
            winter_plant_temp_min = input$winter_plant_temp_min,
            winter_plant_temp_max = input$winter_plant_temp_max
          )
          
          bind_cols(res, meta) %>% mutate(dataset = "Climate")
        }) %>%
        bind_rows()
    }
  })
  
  combined_results <- reactive({
    req(historical_results())
    if (isTruthy(input$climate_file)) {
      bind_rows(historical_results(), future_results())
    } else {
      historical_results()
    }
  })
  
  output$run_parameters_table <- renderTable({
    yrs <- baseline_years_auto()
    data.frame(
      parameter = c("Crop name", "Crop type", "Days to maturity", "Maximum growing days",
                    "Base temperature (°C)", "Optimum temperature (°C)",
                    "Upper temperature cutoff (°C)", "Thermal time method",
                    "Baseline years", "Reference planting date"),
      value = c(input$crop_name, input$crop_type, input$days_to_maturity, input$max_growing_days,
                input$t_base, input$t_opt, input$t_max_cut, input$tt_mode,
                paste0(min(yrs), ":", max(yrs)), input$baseline_planting_mmdd),
      stringsAsFactors = FALSE
    )
  })
  
  output$required_tt_text <- renderPrint({
    req(required_tt_result())
    cat("Crop:", input$crop_name, "\n")
    cat("Crop type:", input$crop_type, "\n")
    cat("Estimated required thermal time:", round(required_tt_result()$required_tt, 2), "\n\n")
    print(required_tt_result()$yearly_required_tt)
  })
  
  output$comparison_summary_table <- renderTable({
    req(combined_results())
    compare_summary(combined_results())
  })
  
  output$results_table <- renderDT({
    req(combined_results())
    datatable(combined_results(), options = list(pageLength = 20, scrollX = TRUE))
  })
  
  
  output$summary_text <- renderText({
    req(combined_results())
    make_summary_text(combined_results())
  })
  
  output$plot_ridge_boxinline <- renderPlot({
    req(combined_results())
    plot_growing_season_ridges(combined_results())
  })
  
  output$plot_monthly_temperature_cycle <- renderPlot({
    req(weather_data())
    if (isTruthy(input$climate_file)) {
      plot_monthly_temperature_cycle(weather_data(), climate_data())
    } else {
      plot_monthly_temperature_cycle(weather_data(), NULL)
    }
  })
  
  output$plot_tmin_tmax_ribbon <- renderPlot({
    req(weather_data())
    if (isTruthy(input$climate_file)) {
      plot_tmin_tmax_ribbon(weather_data(), climate_data())
    } else {
      plot_tmin_tmax_ribbon(weather_data(), NULL)
    }
  })
  
  output$plot_monthly_tt_cycle <- renderPlot({
    req(weather_data())
    if (isTruthy(input$climate_file)) {
      plot_monthly_tt_cycle(weather_data(), climate_data(),
                            t_base = input$t_base, t_opt = input$t_opt,
                            t_max_cut = input$t_max_cut, tt_mode = input$tt_mode)
    } else {
      plot_monthly_tt_cycle(weather_data(), NULL,
                            t_base = input$t_base, t_opt = input$t_opt,
                            t_max_cut = input$t_max_cut, tt_mode = input$tt_mode)
    }
  })
  
  output$plot_annual_tt_boxplot <- renderPlot({
    req(weather_data())
    if (isTruthy(input$climate_file)) {
      plot_annual_tt_boxplot(weather_data(), climate_data(),
                             t_base = input$t_base, t_opt = input$t_opt,
                             t_max_cut = input$t_max_cut, tt_mode = input$tt_mode)
    } else {
      plot_annual_tt_boxplot(weather_data(), NULL,
                             t_base = input$t_base, t_opt = input$t_opt,
                             t_max_cut = input$t_max_cut, tt_mode = input$tt_mode)
    }
  })
  
  output$plot_annual_temperature_boxplot <- renderPlot({
    req(weather_data())
    if (isTruthy(input$climate_file)) {
      plot_annual_temperature_boxplot(weather_data(), climate_data())
    } else {
      plot_annual_temperature_boxplot(weather_data(), NULL)
    }
  })
  
  output$plot_temperature_heatmap <- renderPlot({
    req(weather_data())
    if (isTruthy(input$climate_file)) {
      plot_growing_season_heatmap(weather_data(), climate_data())
    } else {
      plot_growing_season_heatmap(weather_data(), NULL)
    }
  })
  
  output$weather_preview <- renderDT({
    req(weather_data())
    datatable(head(weather_data(), 50), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$weather_info <- renderPrint({
    req(weather_data())
    w <- weather_data()
    cat("Rows:", nrow(w), "\n")
    cat("Period:", as.character(min(w$date)), "to", as.character(max(w$date)), "\n")
  })
  
  output$climate_preview <- renderDT({
    req(input$climate_file)
    datatable(head(climate_data(), 50), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$climate_info <- renderPrint({
    req(input$climate_file)
    w <- climate_data()
    cat("Rows:", nrow(w), "\n")
    cat("Period:", as.character(min(w$date)), "to", as.character(max(w$date)), "\n")
    extra_cols <- setdiff(names(w), c("date", "tmin", "tmax", "year", "doy", "tmean"))
    if (length(extra_cols) > 0) cat("Additional columns:", paste(extra_cols, collapse = ", "), "\n")
  })
}

shinyApp(ui, server)
