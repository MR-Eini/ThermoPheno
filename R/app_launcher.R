#' Launch the ThermoPheno Shiny application
#'
#' Starts the ThermoPheno interactive application for historical analysis and
#' climate change impact assessment of crop phenology.
#'
#' @return The function is called for its side effect of launching the Shiny app.
#' @examples
#' \dontrun{
#' ThermoPheno()
#' }
ThermoPheno <- function() {
  message(
"
========================================
 ThermoPheno
 Thermal-time phenology simulation tool
========================================

Launching interactive app...

Use:
- Upload historical data
- Optionally add climate scenarios
- Explore growing season shifts
")

  app_dir <- system.file("app", package = "ThermoPheno")
  if (app_dir == "") {
    stop("App not found. Please reinstall the package.", call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
