# Winter wheat harvest metrics on the station-years paired in BOTH runs
# (0.1.1 had 8 failed seasons without a simulated harvest date).
#   Rscript validation/compare_common_cases.R
o <- read.csv("pkgdown/assets/validation/archive/v0.1.1/validation_pairs.csv", colClasses = "character")
n <- read.csv("pkgdown/assets/validation/validation_pairs.csv", colClasses = "character")
stopifnot(identical(paste(o$crop_name, o$station_id_pheno, o$year), paste(n$crop_name, n$station_id_pheno, n$year)))
m <- function(obs, sim) {
  e <- as.numeric(as.Date(sim) - as.Date(obs))
  oy <- as.POSIXlt(as.Date(obs))$yday + 1; sy <- as.POSIXlt(as.Date(sim))$yday + 1
  c(n = length(e), bias = mean(e), mae = mean(abs(e)), rmse = sqrt(mean(e^2)), r2 = cor(oy, sy)^2)
}
for (crop in unique(o$crop_name)) {
  k <- o$crop_name == crop & !is.na(o$simulated_harvest_date) & o$simulated_harvest_date != "" &
       !is.na(n$simulated_harvest_date) & n$simulated_harvest_date != ""
  cat(crop, "harvest, common paired station-years\n")
  print(round(rbind(v0.1.1 = m(o$observed_harvest_date[k], o$simulated_harvest_date[k]),
                    v0.2.0 = m(n$observed_harvest_date[k], n$simulated_harvest_date[k])), 3))
}
