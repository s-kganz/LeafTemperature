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

library(tidyverse)
library(foreach)
source("scripts/tower_util.R")

# Read data ----
# Site LAI & metadata
site_meta <- read_csv("data_working/neon_site_metadata.csv") %>%
  left_join(
    read_csv("data_out/neon_sampled_dhp_lai.csv"),
    by=c("site_neon"="site")
  ) %>%
  mutate(
    lai_best = ifelse(is.na(lai), L_dhp_corrected, lai)
  )

# LAD profiles
site_lad <- read_csv("data_out/neon_lad_profiles.csv")

# Flux timestamps
site_flux <- read_csv("data_out/cross_site_flux_partition_qc.csv") %>%
  select(site, timeBgn)

# Set variables to interpolate
regexes <- c("TA_\\d_\\d_\\d", "H2O.*_\\d_\\d_\\d", "CO2_.*\\d_\\d_\\d",
             "WS_\\d_\\d_\\d")
prefixes <- c("TA_", "H2O_", "CO2_", "WS_")

# Interpolation function ----
interpolation_driver <- function(site_amf, site_neon, tower, l) {
  ## Determine z positions corresponding to l ----
  this_lad_profile <- site_lad %>% filter(site == site_neon)
  
  interp_z <- approx(this_lad_profile$cum_lai, this_lad_profile$z, xout=l, rule=2)$y
  
  cat("LAI -> z crosswalk:\n")
  print(data.frame(l=l, z=interp_z))
  
  ## Do interpolation ----
  interp_data <- foreach(pattern=regexes, prefix=prefixes, .combine=cbind) %do% {
    cat("Interpolating regex", pattern, "...\n")
    vars_heights <- get_var_heights(names(tower), pattern, site_amf)
    vars <- vars_heights$vars
    heights <- vars_heights$heights
    interpolate_tower_data(tower, vars, heights, interp_z, name_prefix=prefix)
  }

  stopifnot(nrow(interp_data) == nrow(tower))
  
  ## Convert to long-format table
  interp_long <- interp_data %>%
    mutate(TIMESTAMP = tower$TIMESTAMP,
           site_amf = site_amf,
           site_neon = site_neon) %>%
    pivot_longer(-c(TIMESTAMP, site_amf, site_neon)) %>%
    mutate(varname = str_split_i(name, "_", 1),
           z = str_split_i(name, "_", 2),
           layer_l = l[match(z, interp_z)]) %>%
    select(TIMESTAMP, site_amf, site_neon, layer_l, z, varname, value) %>%
    pivot_wider(names_from=varname, values_from=value)
}

# Run interpolation ----
all_interps <- foreach(i=1:nrow(site_meta), .combine=rbind) %do% {
  site_neon <- site_meta$site_neon[i]
  site_amf  <- site_meta$site_ameriflux[i]
  site_lai  <- site_meta$lai_best[i]
  
  this_l <- seq(0, floor(site_lai)-1)
  
  this_tower <- read_csv(
    file.path("data_working", "neon_flux", str_c(site_amf, ".csv"))
  )
  
  cat("Now processing site", site_neon, "\n")
  # Filter to timestamps that have valid fluxes. This is optional, just makes
  # interpolation go much faster.
  usable_timestamps <- site_flux %>%
    filter(site == site_neon) %>%
    pull(timeBgn)
  
  cat("Usable timestamps:", length(usable_timestamps), "\n")
  
  this_tower_filter <- this_tower %>%
    filter(TIMESTAMP %in% usable_timestamps)
  
  interpolation_driver(site_amf, site_neon, this_tower_filter, this_l)
}

write_csv(all_interps, "data_out/cross_site_interpolated_meteorology.csv")

  