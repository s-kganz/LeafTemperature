#' Interpolate flux tower micrometeorology data at arbitrary z positions.
#'
#' @param tower Data frame containing relevant tower observations.
#' @param vars Vector of names in `tower` to interpolate.
#' @param var_z Vector of heights above ground where `vars` are measured.
#' @param interp_z Vector of heights above ground at which to interpolate `vars`.
#' @param name_prefix String prepended to names of the resulting data frame,
#' which by default are `interp_z`.
#'
#' @return Dataframe with the observations in `tower` interpolated at `interp_z`
#' 
#'
#'
interpolate_tower_data <- function(tower, vars, var_z, interp_z, name_prefix="") {
  stopifnot(length(vars) == length(var_z))
  stopifnot(all(vars %in% names(tower)))
  
  if (!nrow(tower)) {
    warning("Tower data frame is empty! Returning an empty data frame.")
    ret <- data.frame(matrix(nrow=0, ncol=length(interp_z)))
    names(ret) <- str_c(name_prefix, interp_z)
    return(ret)
  }
  
  var_mtx <- tower %>%
    select(all_of(vars)) %>%
    as.matrix()
  
  var_interp <- var_mtx %>%
    apply(
      MARGIN = 1,
      FUN = function(row) {
        tryCatch(
          # This will fire a warning on duplicated x values. This is not a
          # problem so we silence it.
          approx(
            var_z,
            row,
            xout = interp_z,
            rule = 2,
            na.rm = TRUE
          )$y,
          error = function(e)
            return(rep(NA, length(interp_z)))
        )
      }
    )
  
  # If length(interp_z) == 1, var_interp is a vector, otherwise it is a matrix.
  # Transpose only if it is a matrix.
  if (!is.null(dim(var_interp))) {var_interp <- t(var_interp)}
  
  var_interp %>%
    as.data.frame() %>% 
    setNames(str_c(name_prefix, interp_z))
}

# Make a copy of this function in this package so we can memoize it.
# See .onLoad() in zzz.R

#' Return Ameriflux variable info (possibly from cache). This function only exists so that it can be memoized.
#'
#' @return see \link[amerifluxr]{amf_var_info}
#' @keywords internal
#' @importFrom memoise memoise
#' @importFrom amerifluxr amf_var_info
amf_var_heights <- function() {amerifluxr::amf_var_info()}

#' Get variables matching a regex and their height of observation.
#'
#' @param all_vars Character vector of variables in a tower data frame.  
#' @param regex Character string used to `grep` variables out of `all_vars`.
#' @param site Site code used to filter results from `amf_var_info`.
#'
#' @return a list containing the regex'd varnames and their heights
#' @keywords internal
#'
#'
get_var_heights <- function(all_vars, regex, site) {
  if (!exists("amf_var_heights_")) {amf_var_heights_ <- amf_var_heights()}
  stopifnot(site %in% amf_var_heights_$Site_ID)
  
  vars_select <- grep(regex, all_vars, value=TRUE)
  if (!length(vars_select)) warning("No matching variables found in table!")
  
  amf_var_heights_site <- amf_var_heights_ %>% filter(Site_ID == site)
  if (!nrow(amf_var_heights_site)) warning("Variable heights not found!")
  
  vars_height <- amf_var_heights_site$Height[
    match(vars_select, amf_var_heights_site$Variable)
  ]
  
  list(
    vars=vars_select,
    heights=vars_height
  )
} 

write_if_not_exist <- function(table, path, ...) {
  if (!file.exists(path)) {
    write_csv(table, path, ...)
  } else {
    warning(paste("File", path, "already exists. Table not written."))
  }
}

if (sys.nframe() == 0) {
  a_tower <- read_csv("data_working/neon_flux/US-xWR.csv")
  
  interp_z <- c(0, 6, 10, 15)
  
  # Test interpolation of [CO2], [H2O], Ta, wind speed
  patterns <- c("TA_\\d_\\d_\\d", "WS_\\d_\\d_\\d", "CO2_.*\\d_\\d_\\d",
                "H2O_.*\\d_\\d_\\d")
  prefixes <- c("TA", "WS", "CO2", "H2O")
  interp_data <- foreach(pattern=patterns, prefix=prefixes, .combine=cbind) %do% {
    vars_heights <- get_var_heights(names(a_tower), pattern, "US-xWR")
    vars <- vars_heights$vars
    heights <- vars_heights$heights
    
    interpolate_tower_data(a_tower, vars, heights, interp_z, name_prefix=prefix)
  }
}
