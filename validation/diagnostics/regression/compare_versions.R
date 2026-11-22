# Compare saved outputs of run_exported_functions.R for two versions.
#   Rscript validation/diagnostics/regression/compare_versions.R old.rds new.rds
a <- commandArgs(trailingOnly = TRUE)
old <- readRDS(a[[1]]); new <- readRDS(a[[2]])
cat("Old version:", old$version, " New version:", new$version, "\n\n")
keys <- setdiff(names(old), "version")
rows <- list()
for (k in keys) for (f in names(old[[k]])) {
  o <- old[[k]][[f]]; n <- new[[k]][[f]]
  if (f == "run_simulation") {
    o_cmp <- o; n_cmp <- n
    dates <- c("planting_date", "maturity_date", "harvest_date")
    same_dates <- all(vapply(dates, function(d) identical(as.character(o[[d]]), as.character(n[[d]])), logical(1)))
    same_status <- identical(o$status, n$status)
    maxd <- max(abs(o$accumulated_tt - n$accumulated_tt), na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(case = k, fun = f, identical = identical(o, n),
      max_abs_numeric_diff = maxd, same_dates_and_status = same_dates && same_status,
      old_value = paste0("n=", nrow(o), "; mature=", sum(o$status == "mature")),
      new_value = paste0("n=", nrow(n), "; mature=", sum(n$status == "mature")))
  } else {
    ov <- if (is.list(o)) o$required_tt else o
    nv <- if (is.list(n)) n$required_tt else n
    yo <- if (is.list(o)) o$yearly_required_tt[[ncol(o$yearly_required_tt)]] else NA
    yn <- if (is.list(n)) n$yearly_required_tt[[ncol(n$yearly_required_tt)]] else NA
    rows[[length(rows) + 1]] <- data.frame(case = k, fun = f,
      identical = identical(ov, nv) && identical(yo, yn),
      max_abs_numeric_diff = max(c(abs(ov - nv), abs(yo - yn)), na.rm = TRUE),
      same_dates_and_status = NA,
      old_value = sprintf("%.10f", ov), new_value = sprintf("%.10f", nv))
  }
}
tab <- do.call(rbind, rows)
options(width = 220)
print(tab, row.names = FALSE)
