.onLoad <- function(libname, pkgname) {
  # Make %do% an operator for infix syntax
  `%do%` <<- foreach::`%do%`
}
