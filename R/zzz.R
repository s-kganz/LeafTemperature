.onLoad <- function(libname, pkgname) {
  amf_var_heights <<- memoise::memoise(
    # Don't ping ameriflux API for a day
    amf_var_heights, cache=cachem::cache_mem(max_age=86400)
  )
}