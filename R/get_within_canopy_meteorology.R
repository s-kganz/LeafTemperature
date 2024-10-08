# To run the energy balance model, we have to calculate canopy micrometeorology
# at each of the layers. But, each layer is defined by L, *not* height above
# ground. So, we do the following for each site:
# 
# - Set cumulative L's
# - Use data_out/neon_lad_profiles.csv to map L's to heights above ground
# - Interpolate CO2, H2O, TA, and WS at these heights
# 
# Note that this step is quite data-intensive, so we filter each tower dataset
# to timestamps that have valid fluxes.

# Read data ----
# Site LAI & metadata


# Interpolation functions ----
neon_interpolation_driver <- function(site_amf, site_neon, site_lad, tower, l) {
  ## Determine z positions corresponding to l ----
  this_lad_profile <- site_lad %>% filter(site == site_neon)
  
  interp_z <- approx(this_lad_profile$cum_lai, this_lad_profile$z, xout=l, rule=2)$y
  
  message(paste("Interpolating", site_neon))
  message("LAI -> z crosswalk:")
  message(paste0(capture.output(data.frame(l=l, z=interp_z)), collapse = "\n"))

  ## Do interpolation ----
  # Set variables to interpolate
  regexes <- c("TA_\\d_\\d_\\d", "H2O.*_\\d_\\d_\\d", "CO2_.*\\d_\\d_\\d",
               "WS_\\d_\\d_\\d")
  
  prefixes <- c("TA_", "H2O_", "CO2_", "WS_")
  interp_data <- foreach(pattern=regexes, prefix=prefixes, .combine=cbind) %do% {
    message(paste("Interpolating regex", pattern))
    vars_heights <- get_var_heights(names(tower), pattern, site_amf)
    vars <- vars_heights$vars
    heights <- vars_heights$heights
    suppressWarnings(interpolate_tower_data(tower, vars, heights, interp_z, name_prefix=prefix))
  }

  stopifnot(nrow(interp_data) == nrow(tower))
  
  ## Convert to long-format table
  interp_long <- interp_data %>%
    mutate(TIMESTAMP = tower$TIMESTAMP,
           site_amf = site_amf,
           site_neon = site_neon) %>%
    pivot_longer(-c(TIMESTAMP, site_amf, site_neon)) %>%
    mutate(varname = stringr::str_split_i(name, "_", 1),
           z = stringr::str_split_i(name, "_", 2),
           layer_l = l[match(z, interp_z)]) %>%
    select(TIMESTAMP, site_amf, site_neon, layer_l, z, varname, value) %>%
    pivot_wider(names_from=varname, values_from=value)
}

do_neon_interpolation <- function(site_meta, site_lai, site_flux, site_lad, outdir) {
  # Run interpolation
  all_interps <- foreach(i=1:nrow(site_meta), .combine=rbind) %do% {
    site_neon <- site_meta$site_neon[i]
    site_amf  <- site_meta$site_ameriflux[i]
    site_lai  <- site_meta$lai_best[i]
    
    this_l <- seq(0, floor(site_lai))
    
    this_tower <- read_csv(
      file.path(outdir, "amf_tower_data", str_c(site_amf, ".csv")),
      col_types = cols()
    )
    
    #cat("Now processing site", site_neon, "\n")
    # Filter to timestamps that have valid fluxes. This is optional, just makes
    # interpolation go much faster.
    usable_timestamps <- site_flux %>%
      filter(site == site_neon) %>%
      pull(timeBgn)
    
    #cat("Usable timestamps:", length(usable_timestamps), "\n")
    
    this_tower_filter <- this_tower %>%
      filter(TIMESTAMP %in% usable_timestamps)
    
    neon_interpolation_driver(site_amf, site_neon, site_lad, this_tower_filter, this_l)
  }
  
  write_if_not_exist(all_interps, 
                     file.path(outdir, "cross_site_interpolated_meteorology.csv"))
}

do_old_wref_interpolation <- function() {
  # Interpolation in the spreadsheet has to be manually specified
  wref_flux <- read_csv("data_out/old_wref_clean.csv")
  
  ta_vars <- c("TA_02", "TA_70")
  rh_vars <- c("RH_02", "RH_70")
  heights <- c(2, 70)
  
  # Get heights at which to interpolate
  lai <- site_meta$lai_best[match("WREF", site_meta$site_neon)]
  interp_l <- seq(0, floor(lai))
  
  this_lad_profile <- site_lad %>% filter(site == "WREF")
  
  interp_z <- approx(this_lad_profile$cum_lai, this_lad_profile$z, 
                     xout=interp_l, rule=2)$y
  
  interp_ta <- interpolate_tower_data(
    wref_flux, ta_vars, heights, interp_z, name_prefix="TA_"
  )
  interp_rh <- interpolate_tower_data(
    wref_flux, rh_vars, heights, interp_z, name_prefix="RH_"
  )
  
  wref_interp <- cbind(interp_ta, interp_rh) %>%
    mutate(TIMESTAMP=wref_flux$TIMESTAMP) %>%
    pivot_longer(-c(TIMESTAMP)) %>%
    mutate(varname = stringr::str_split_i(name, "_", 1),
           z = stringr::str_split_i(name, "_", 2),
           layer_l = interp_l[match(z, interp_z)]) %>%
    select(TIMESTAMP, layer_l, z, varname, value) %>%
    pivot_wider(names_from=varname, values_from=value)
  
  write_if_not_exist(wref_interp, "data_out/old_wref_interpolated_meteorology.csv")
}

#' Interpolate within-canopy meteorology to model layers
#'
#' @param site_meta Site metadata table \code{\link{site_meta}}
#' @param site_lai Site leaf area index table \code{\link{modis_lai}}
#' @param site_lad Site leaf area density profiles generated by [fit_neon_lidar_k()].
#' @param site_flux Partitioned eddy covariance flux generated by [partition_neon_flux()].
#' @param outdir Output directory where interpolated meteorology table will be written.
#'
#' @description
#' In all sites, model layers are not necessarily at the same physical height as meteorological sensors. This is because model layers are based on cumulative LAI, which does not change linearly with height. Therefore, we use linear interpolation to approximate microclimate conditions at the position of each model layer. This script looks up the sensor heights of relevant variabels in the tower data and then interpolates to the layer positions.
#' 
#'
#' @export
#'
get_within_canopy_meteorology <- function(site_meta, site_lai, site_lad, site_flux, outdir) {
  site_meta <- site_meta %>%
    left_join(
      site_lai,
      by=c("site_ameriflux"="ID")
    ) %>%
    mutate(
      lai_best = ifelse(is.na(lai), mod_lai, lai)
    )
  
  # Only pull the timestamps from the flux table so we cut down on processing
  # time.
  site_flux <- site_flux %>%
    select(site, timeBgn)
  
  do_neon_interpolation(site_meta, site_lai, site_flux, site_lad, outdir)
}
  
