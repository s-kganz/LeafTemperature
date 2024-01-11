# This file joins tower micrometeorology with T_canopy QC and uncertainty
# estimates from NEON

# Libraries ----

library(tidyverse)
library(bigleaf)
library(lubridate)
library(amerifluxr)
library(neonUtilities)
source("scripts/energy_balance/leaf_eb_still.R")

var_info  <- amf_var_info()
site_info <- read_csv("data_working/neon_site_metadata.csv")

# Constants ----

# Simple top-of-canopy measurements that do not require interpolation
toc_vars <- c("TIMESTAMP", "P", "PA", "LE", "SW_DIF", "SW_OUT", "LW_IN", 
              "LW_OUT", "USTAR")

# Regexes to collect variables in replicate or in a vertical profile
sw_in_regex <- "SW_IN_\\d_1_\\d"
ta_regex    <- "TA_\\d_\\d_\\d"
# Horizontal qualifier 1 distinguishes tower sensors from those in soil plots
tc_regex    <- "T_CANOPY_1_\\d_\\d"
ws_regex    <- "WS_\\d_\\d_\\d"
mr_regex    <- "H2O_MIXING_RATIO_\\d_\\d_\\d"
g_regex     <- "G_\\d_\\d_\\d"
ts_regex    <- "TS_\\d_1_\\d"
swc_regex   <- "SWC_\\d_\\d_\\d"
par_regex   <- "PPFD_IN_\\d_\\d_\\d"
par_toc_regex <- "PPFD_IN_\\d_1_\\d"

# Clean and interpolate tower micrometeorology ----
# Canopy temperature and micrometorology data are generally not measured
# at the same point. So we use interpolation to get data for the same location
# in space.

get_site_tc_data <- function(site, utc_offset) {
  data <- read_csv(file.path("data_working", "neon_flux", str_c(site, ".csv")),
                   col_types=cols())
  
  cat(site, "\n")
  cat(nrow(data), "initial observations", "\n")
  
  sw_in_vars <- grep(sw_in_regex, names(data), value=TRUE)
  ta_vars    <- grep(ta_regex, names(data), value=TRUE)
  tc_vars    <- grep(tc_regex, names(data), value=TRUE)
  ws_vars    <- grep(ws_regex, names(data), value=TRUE)
  mr_vars    <- grep(mr_regex, names(data), value=TRUE)
  g_vars     <- grep(g_regex, names(data), value=TRUE)
  ts_vars    <- grep(ts_regex, names(data), value=TRUE)
  swc_vars   <- grep(swc_regex, names(data), value=TRUE)
  par_vars   <- grep(par_regex, names(data), value=TRUE)
  par_toc_vars <- grep(par_toc_regex, names(data), value=TRUE)
  
  cat("TOC SW_IN:", str_c(sw_in_vars, sep=" "), "\n")
  cat("TA:", str_c(ta_vars, sep=" "), "\n")
  cat("TC:", str_c(tc_vars, sep=" "), "\n")
  cat("WS:", str_c(ws_vars, sep=" "), "\n")
  cat("MIX:", str_c(mr_vars, sep=" "), "\n")
  cat("G", str_c(g_vars, sep=" "), "\n")
  cat("Surface TS:", str_c(ts_vars, sep=" "), "\n")
  cat("Surface SWC:", str_c(swc_vars, sep=" "), "\n")
  cat("PAR:", str_c(par_vars, sep=" "), "\n")
  cat("TOC PAR:", str_c(par_toc_vars, sep=" "), "\n")
  
  # Determine heights for interpolated vars ahead of time so we can interpolate
  # without having to join anything (which is slow).
  this_heights <- var_info %>% filter(Site_ID == site)
  ta_heights   <- this_heights$Height[match(ta_vars, this_heights$Variable)]
  ws_heights   <- this_heights$Height[match(ws_vars, this_heights$Variable)]
  tc_heights   <- this_heights$Height[match(tc_vars, this_heights$Variable)]
  par_heights  <- this_heights$Height[match(par_vars, this_heights$Variable)]
  
  # Some sites have not updated the var info table with the fact that H2O conc is
  # now reported as mixing ratios. If the result is NA, try again with 
  # the old variable names.
  mr_heights <- this_heights$Height[match(mr_vars, this_heights$Variable)]
  
  if (all(is.na(mr_heights))) {
    cat("No heights found for mixing ratio, trying again with old var name.", "\n")
    mr_vars <- str_replace(mr_vars, "H2O_MIXING_RATIO", "H2O")
    mr_heights <- this_heights$Height[match(mr_vars, this_heights$Variable)]
  }
  
  if (!all(is.na(mr_heights))) {
    cat("Found mixing ratio heights with old names!", "\n")
  } else {
    cat("No mixing heights with old names either!", "\n")
  }
  
  # US-xTE has not updated the var metadata table with measurement heights for
  # T_CANOPY nor H2O conc. These values are taken from the NEON releases
  # of related data.
  if (site == "US-xTE") {
    tc_heights <- c(
      26.86, 18.64, 10.11, 3.42, 0.9
    )
    mr_heights <- c(
      59.17, 59.17, 59.17, 32.71, 32.71, 26.85, 26.85, 18.64, 18.64,
      10.11, 10.11, 3.42, 3.42, 0.90, 0.90
    )
  }
  
  data_select <- data %>%
    select(all_of(toc_vars), matches(ta_regex), matches(tc_regex),
           matches(mr_regex), matches(ws_regex), matches(par_regex)) %>%
    mutate(
      SW_IN_TOC = rowMeans(data %>% select(matches(sw_in_regex)), na.rm=TRUE),
      PPFD_IN_TOC = rowMeans(data %>% select(matches(par_toc_regex)), na.rm=TRUE),
      G = rowMeans(data %>% select(matches(g_regex)), na.rm=TRUE),
      TS = rowMeans(data %>% select(matches(ts_regex)), na.rm=TRUE),
      SWC = rowMeans(data %>% select(matches(swc_regex)), na.rm=TRUE)
    ) %>%
    nest(
      tc = matches(tc_regex),
      ta = matches(ta_regex),
      ws = matches(ws_regex),
      mr = matches(mr_regex),
      par = matches(par_regex)
    )
  
  data_interp <- data_select %>%
    # Drop observations that don't have at least two non-NAs for interpolation
    filter(
      map_int(tc, ~sum(!is.na(.x))) >= 2,
      map_int(ta, ~sum(!is.na(.x))) >= 2,
      map_int(ws, ~sum(!is.na(.x))) >= 2,
      map_int(mr, ~sum(!is.na(.x))) >= 2,
      map_int(par, ~sum(!is.na(.x))) >= 2
    ) %>%
    mutate(
      tc = pmap(
        list(tc, ta, ws, mr, par),
        function(tc, ta, ws, mr, par) {
          tc %>%
            pivot_longer(everything(), values_to = "t_canopy") %>%
            mutate(
              z = tc_heights,
              zmax = max(z),
              ta_interp = spline(ta_heights, unlist(ta), xout=z)$y,
              ws_interp = spline(ws_heights, unlist(ws), xout=z)$y,
              par_interp = spline(par_heights, unlist(par), xout=z)$y,
              # Mixing ratio instruments have duplicated z coordinates, spline
              # fires a warning for this but it really doesn't matter.
              mr_interp = suppressWarnings(spline(mr_heights, unlist(mr), xout=z)$y)
            )
        }
      )
    ) %>%
    select(-ta, -ws, -mr, -par) %>%
    unnest(tc) %>%
    # Model blows up if wind == 0 so drop those observations
    filter(ws_interp >= 0)
    
  cat(nrow(data_interp), "TC observations before filtering", "\n")
  
  data_final <- data_interp %>%
    # Compute relative humidity
    mutate(
      # Mixing ratio is provided by ameriflux in mmol H2O/mol dry air so divide
      # by 1000 to convert to mol ratio
      e    = PA * (1 + (mr_interp/1000)^-1)^-1,
      esat = saturation_vapor_pressure(ta_interp + 273.15),
      rh   = e / esat,
      vpd  = esat - e,
      # Clean up the G column
      G = replace_na(G, 0)
    ) %>%
    # Filter to reasonable values of RH and VPD
    filter(rh >= 0, rh <= 1, vpd < 6) %>%
    # Growing season only
    filter(month(TIMESTAMP) %in% c(5:9)) %>%
    # Column cleanup
    mutate(
      SITE = site
    ) %>%
    select(TIMESTAMP, SITE, everything(), -name)
  
  cat(nrow(data_final), "TC observations after filtering", "\n")
  cat("\n\n")
  
  data_final
}

all_amf_data <- bind_rows(lapply(
  1:nrow(site_info), function(i) {
    site_amf <- site_info$site_ameriflux[i]
    site_neon <- site_info$site_neon[i]
    utc_offset <- site_info$utc_offset[i]
    get_site_tc_data(site_amf, utc_offset) %>% mutate(SITE_NEON=site_neon)
  }
))

# Join NEON QC and uncertainty ----

# Now pull NEON IRBT data to attach QC flag and uncertainty. Join schema
# is on site ID, timestamp, and z position. We keep the tc value from the
# Ameriflux portal to verify that the join is correct.
token <- readChar("data_in/neon_token", file.info("data_in/neon_token")$size)

prod <- loadByProduct(
  "DP1.00005.001",
  site=site_info$site_neon,
  tabl="IRBT_30_minute",
  token=token,
  nCores=4,
  check.size=FALSE
)

irbt <- prod$IRBT_30_minute

irbt_data <- prod$IRBT_30_minute %>%
  select(siteID, horizontalPosition, verticalPosition,
         startDateTime, bioTempMean, bioTempExpUncert,
         finalQF) %>%
  filter(finalQF == 0) %>%
  mutate(HOR.VER = str_c(horizontalPosition, verticalPosition, sep="."))

# Some vertical positions are duplicated and the zoffsets vary by < 1m. Average
# these together to silence later warnings.
irbt_pos <- prod$sensor_positions_00005 %>%
  select(siteID, HOR.VER, zOffset)

irbt_join <- irbt_data %>%
  left_join(irbt_pos,
            by=c("siteID", "HOR.VER"))

amf_neon_joined <- all_amf_data %>%
  inner_join(irbt_join,
             by=c("SITE_NEON"="siteID", "TIMESTAMP"="startDateTime",
                  "z"="zOffset")) %>%
  filter(finalQF == 0) %>%
  select(-HOR.VER,
         -finalQF) %>%
  rename(t_canopy_exp_uncertainty = "bioTempExpUncert",
         t_canopy_neon = "bioTempMean") %>%
  select(TIMESTAMP, SITE, SITE_NEON, t_canopy, t_canopy_neon, 
         t_canopy_exp_uncertainty, everything())

stopifnot(all.equal(amf_neon_joined$t_canopy, amf_neon_joined$t_canopy_neon))

write_csv(amf_neon_joined, "data_out/cross_site_tc_data.csv")
