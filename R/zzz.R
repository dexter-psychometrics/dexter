
# this avoids having to attach shinyBS

.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))
}

.onAttach <- function(libname, pkgname) {
  shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))
}

.onUnload <- function (libpath) {
  library.dynam.unload("dexter", libpath)
}