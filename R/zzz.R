.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
"ThermoPheno loaded.

Run the app using:
  ThermoPheno()

Example files are available via system.file('extdata', package = 'ThermoPheno')."
  )
}
