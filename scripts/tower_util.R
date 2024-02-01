library(tidyverse)
library(amerifluxr)

if (!exists("amf_var_heights_")) {amf_var_heights_ <- amf_var_info()}

#' Interpolate flux tower micrometeorology data at arbitrary z positions.
#'
#' @param tower Data frame containing relevant tower observations.
#' @param vars Vector of names in `tower` to interpolate.
#' @param var_z Vector of heights above ground where `vars` are measured.
#' @param interp_z Vector of heights above ground at which to interpolate `vars`.
#' @param name_prefix String prepended to names oft he resulting data frame,
#' which by default are `interp_z`.
#'
#' @return
#' @export
#'
#' @examples
interpolate_tower_data <- function(tower, vars, var_z, interp_z, name_prefix="") {
  stopifnot(length(vars) == length(var_z))
  stopifnot(all(vars %in% names(tower)))
  
  var_mtx <- tower %>%
    select(all_of(vars)) %>%
    as.matrix()
  
  var_interp <- var_mtx %>%
    apply(MARGIN=1, 
          FUN=function(row) {
            tryCatch(
              # This will fire a warning on duplicated x values. This is not a
              # problem so we silence it.
              approx(var_z, row, xout=interp_z, rule=2, na.rm=TRUE)$y, 
              error=function(e) return(rep(NA, length(interp_z))))
          }
    ) %>% t()
  
  var_interp %>%
    as.data.frame() %>% 
    setNames(str_c(name_prefix, interp_z))
}

#' Get variables matching a regex and their height of observation.
#'
#' @param all_vars Character vector of variables in a tower data frame.  
#' @param regex Character string used to `grep` variables out of `all_vars`.
#' @param site Site code used to filter results from `amf_var_info`.
#'
#' @return
#' @export
#'
#' @examples
get_var_heights <- function(all_vars, regex, site) {
  if (!exists("amf_var_heights")) {amf_var_heights_ <- amf_var_info()}
  stopifnot(site %in% amf_var_heights_$Site_ID)
  
  vars_select <- grep(regex, all_vars, value=TRUE)
  
  amf_var_heights_site <- amf_var_heights_ %>% filter(Site_ID == site)
  
  vars_height <- amf_var_heights_site$Height[
    match(vars_select, amf_var_heights_site$Variable)
  ]
  
  list(
    vars=vars_select,
    heights=vars_height
  )
} 

if (sys.nframe() == 0) {
  a_tower <- read_csv("data_working/neon_flux/US-xWR.csv")
  
  interp_z <- c(0, 6, 10, 15)
  
  # Test interpolation of [CO2], [H2O], Ta, wind speed
  patterns <- c("TA_\\d_\\d_\\d", "WS_\\d_\\d_\\d", "CO2_.*\\d_\\d_\\d",
                "H2O_.*\\d_\\d_\\d")
  prefixes <- c("TA", "WS", "CO2", "H2O")
  library(foreach)
  interp_data <- foreach(pattern=patterns, prefix=prefixes, .combine=cbind) %do% {
    vars_heights <- get_var_heights(names(a_tower), pattern, "US-xWR")
    vars <- vars_heights$vars
    heights <- vars_heights$heights
    
    interpolate_tower_data(a_tower, vars, heights, interp_z, name_prefix=prefix)
  }
}
