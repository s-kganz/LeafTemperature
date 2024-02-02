# This script runs the final energy balance model.

library(tidyverse)
library(suncalc)
library(photosynthesis)
source("scripts/energy_balance/leaf_eb_numeric.R")
# Read source data and constants ----

## Whole-canopy fluxes ----
flux <- read_csv("data_out/cross_site_flux_partition_qc.csv") %>%
  select(site, timeBgn, GPP_U50_f, data.fluxH2o.nsae.flux) %>%
  rename(SITE_NEON=site,
         TIMESTAMP=timeBgn,
         TOC_GPP=GPP_U50_f,
         TOC_LE=data.fluxH2o.nsae.flux) %>%
  mutate(
    # This process involves rearranging a lot of data. Later we will be grouping
    # model output by timestamp and site. Grouping by a single serial number
    # might be more efficient.
    FLUX_OBS_NUMBER = 1:nrow(.)
  )

## Site Lat, Lon LAI ----
# If an LAI value exists in the literature, use it. Otherwise
# use the corrected LAI derived from DHPs.
site_meta <- read_csv("data_working/neon_site_metadata.csv") %>%
  left_join(
    read_csv("data_out/neon_sampled_dhp_lai.csv"),
    by=c("site_neon"="site")
  ) %>%
  mutate(
    lai_best = ifelse(is.na(lai), L_dhp_corrected, lai)
  ) %>%
  select(site_neon, site_ameriflux, tower_lat, tower_lon, lai_best) %>%
  rename(
    SITE_NEON=site_neon,
    SITE_AMF=site_ameriflux,
    TOWER_LAT=tower_lat,
    TOWER_LON=tower_lon,
    LAI=lai_best
  ) %>%
  # Calculate layer positions for each site based on LAI. We do this in layers
  # with LAI=1 because this makes the total canopy gas exchange equal to the 
  # simple sum of all the layers.
  mutate(
    LAYER_L = map(LAI, function(x) seq(0, floor(x)-1))
  )

## Photosynthesis submodel constants ----
# Add variable prefixes to reduce chance of confusion.
aq_constants <- read_csv("data_out/cross_site_aq_constants.csv") %>%
  rename_with(
    .fn=function(x) str_c("AQ_", x),
    .cols=c(-site)
  ) %>%
  rename(SITE_NEON=site) %>%
  select(SITE_NEON, everything())

medlyn_constants <- read_csv("data_out/cross_site_medlyn_coefficients.csv") %>%
  rename_with(
    .fn=function(x) str_c("MEDLYN_", x),
    .cols=c(-site)
  ) %>%
  rename(SITE_AMF=site)

## Radiation submodel constants ----
lidar_constants <- read_csv("data_out/neon_lidar_constants.csv") %>%
  rename_with(
    .fn=function(x) str_c("LIDAR_", x),
    .cols=c(-site)
  ) %>%
  rename(SITE_NEON=site,
         LIDAR_k_lai=LIDAR_lai_coef) %>%
  select(SITE_NEON, everything())

## Top of canopy radiation fluxes ----
tower_files <- list.files("data_working/neon_flux/", pattern="*.csv",
                          full.names=TRUE)

tower_toc_rad <- lapply(tower_files, function(file) {
  site_amf <- basename(file) %>% str_remove(".csv")
  
  this_tower <- read_csv(file, col_types=cols())
  
  this_tower %>%
    select(TIMESTAMP, PA, SW_DIF, SW_OUT, LW_IN, LW_OUT) %>%
    # Indicate that the fluxes are top-of-canopy
    rename_with(
      .fn=function(x) str_c("TOC_", x),
      .cols=-TIMESTAMP
    ) %>%
    mutate(
      TOC_SW_IN = rowMeans(this_tower %>% select(matches("SW_IN_1_1_\\d")), 
                            na.rm=TRUE),
      TOC_TA = rowMeans(this_tower %>% select(matches("TA_1_1_\\d")),
                        na.rm=TRUE),
      SITE_AMF=site_amf
    ) %>%
    # Model isn't meaningful at night
    filter(TOC_SW_IN > 100) %>%
    select(SITE_AMF, TIMESTAMP, everything()) %>%
    # All the fluxes + air temperature need to be present
    drop_na() %>%
    rename_with(
      .fn=function(x) str_c("RAD_", x),
      .cols=-c(TIMESTAMP, SITE_AMF)
    )
}) %>% bind_rows()

## Within-canopy meteorology ----
interp_meteo <- read_csv("data_out/cross_site_interpolated_meteorology.csv") %>%
  rename(
    SITE_AMF=site_amf,
    SITE_NEON=site_neon,
    LAYER_L=layer_l
  ) %>%
  rename_with(
    .fn=function(x) str_c("LAYER_", x),
    .cols=-c(TIMESTAMP, SITE_AMF, SITE_NEON, LAYER_L)
  )

# Join constants with flux table ----
flux_constants <- flux %>%
  inner_join(site_meta, by="SITE_NEON") %>%
  inner_join(aq_constants, by="SITE_NEON") %>%
  inner_join(medlyn_constants, by="SITE_AMF") %>%
  inner_join(lidar_constants, by="SITE_NEON")

# Join top-of-canopy radiation with flux table ----
flux_rad <- flux_constants %>%
  inner_join(tower_toc_rad,
             by=c("SITE_AMF", "TIMESTAMP"))

# Calculate solar position and black-leaf radiation constants ----
flux_solar <- flux_rad %>%
  bind_cols(
    getSunlightPosition(
      data = data.frame(
        date=.$TIMESTAMP,
        lat=.$TOWER_LAT,
        lon=.$TOWER_LON
      )
    )
  ) %>%
  # Cleanup
  select(-date, -lat, -lon) %>%
  rename(LIDAR_sun_altitude=altitude,
         LIDAR_sun_azimuth=azimuth) %>%
  # Calculate black-leaf radiation constants. This assumes spherical leaf angle
  # distribution. See Appendix in Wang and Leuning (1998) and Ch. 6 in Goudriaan
  # and van Laar (1994).
  mutate(
    LIDAR_sun_zenith = (pi/2) - LIDAR_sun_altitude,
    LIDAR_kb_black = -0.5 / cos(LIDAR_sun_zenith),
    LIDAR_kd_black = -0.8 # constant
  )

# Expand layer L and join meteorology ----
flux_layer_l <- flux_solar %>%
  unnest(LAYER_L) %>%
  relocate(LAYER_L, .after=last_col()) %>%
  inner_join(
    interp_meteo,
    by=c("TIMESTAMP"="TIMESTAMP", "SITE_AMF", "SITE_NEON", "LAYER_L")
  )

# TODO should we calculate canopy reflectance coefficients? Below figure
# suggests that 0.1 is a reasonable approximation and that this is unaffected
# by the proportion of diffuse light hitting the canopy.
# flux_layer_l %>%
#   filter(RAD_TOC_SW_IN > 100) %>%
#   mutate(proportion_diffuse = RAD_TOC_SW_DIF / RAD_TOC_SW_IN,
#          proportion_diffuse_bin = round(proportion_diffuse, digits=1),
#          proportion_diffuse_bin = pmax(0, pmin(proportion_diffuse_bin, 1)),
#          canopy_reflectance = RAD_TOC_SW_OUT / RAD_TOC_SW_IN,
#          canopy_reflectance = pmax(0, pmin(canopy_reflectance, 1))) %>%
#   ggplot(aes(x=factor(proportion_diffuse_bin), y=canopy_reflectance)) +
#   geom_boxplot() +
#   scale_y_log10() +
#   facet_wrap(~ SITE_NEON, scales="free_y") +
#   theme_bw()

# Run radiation submodel ----
flux_rad_model <- flux_layer_l %>%
  mutate(
    RAD_SKY_EMIS = RAD_TOC_LW_IN / (eb_constants_$sb * (RAD_TOC_TA+273.15)^4),
    RAD_SKY_EMIS = pmax(0, pmin(1, RAD_SKY_EMIS)),
    RAD_SW_ABS = absorbed_shortwave_radiation(
      SW_beam=pmax(0, RAD_TOC_SW_IN - RAD_TOC_SW_DIF),
      SW_dif=RAD_TOC_SW_DIF,
      l=LAYER_L,
      kb_real=-LIDAR_kb_real,
      kd_real=-LIDAR_kd_real,
      kb=-LIDAR_kb_black,
      kd=-LIDAR_kd_black,
      s=0.5,
      rho_b=0.1,
      rho_d=0.1,
      aslist=FALSE
    ),
    RAD_LW_ABS = absorbed_downwelling_longwave_radiation(
      Ta=(RAD_TOC_TA+273.15),
      l=LAYER_L,
      L=LAI,
      kb=-LIDAR_kb_black,
      ea=RAD_SKY_EMIS,
      ef=0.95,
      es=0.9,
      kd=-LIDAR_kd_black
    )
  )

# Run photosynthesis submodel ----
flux_ps_model <- flux_rad_model %>%
  # Calculate layer contribution to GPP
  mutate(
    PS_LAYER_GPP_WEIGHT = marshall_biscoe_1980(
      RAD_SW_ABS, AQ_k_sat, AQ_phi_J, AQ_theta_J
    )
  ) %>%
  # Sum the contributions within the whole canopy
  group_by(FLUX_OBS_NUMBER) %>%
  mutate(
    PS_SUM_GPP_WEIGHT = sum(PS_LAYER_GPP_WEIGHT)
  ) %>%
  ungroup() %>%
  # Allocate GPP to each layer such that the sum equals the canopy total
  mutate(
    PS_LAYER_GPP = TOC_GPP * (PS_LAYER_GPP_WEIGHT / PS_SUM_GPP_WEIGHT)
  ) %>%
  # Do some calculations on air moisture
  mutate(
    LAYER_SAT_VP = saturation_vapor_pressure(LAYER_TA+273.15),
    # Convert mmol / mol to mol / mol
    LAYER_VP = RAD_TOC_PA * (LAYER_H2O / 1000),
    # Force reasonable values
    LAYER_VPD = pmax(LAYER_SAT_VP - LAYER_VP, 0.1),
    LAYER_RH = 100 * (1 - (LAYER_VPD / LAYER_SAT_VP)),
  ) %>%
  # Drop observations with missing data for getting stomatal conductance
  drop_na(c(MEDLYN_g1, LAYER_CO2, LAYER_H2O)) %>%
  # Calculate stomatal conductance in the layer
  mutate(PS_LAYER_GS = medlyn_gs(PS_LAYER_GPP, LAYER_VPD, LAYER_CO2, MEDLYN_g1, MEDLYN_g0))

# Run final energy balance model ----
flux_eb_result <- flux_ps_model %>%
  mutate(
    EB_MODEL = pmap(
      list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, RAD_TOC_PA, LAYER_RH,
           RAD_SW_ABS, RAD_LW_ABS),
      energy_balance_driver,
      bounds=20,
      run_no_transp=TRUE,
      .progress=TRUE
    )
  ) %>%
  unnest_wider(EB_MODEL, names_sep="_")

#write_csv(flux_eb_result, "data_out/cross_site_eb_result.csv")

#flux_eb_unnest <- read_csv("data_out/cross_site_eb_result.csv")

# Plots! ----
# Tl by site
(plot_tl_hist <- flux_eb_result %>%
   ggplot(aes(x=LAYER_TA, y=EB_MODEL_Tl-273.15-LAYER_TA)) +
   geom_point(aes(color=LAYER_L/LAI), alpha=0.1) +
   #geom_abline(slope=1, intercept=0, color="red") +
   facet_wrap(~ SITE_NEON) +
   scale_color_viridis_c() +
   theme_bw()
)

# Layer-by-layer plot of the two fluxes
(plot_rad_fluxes <- flux_ps_model %>%
  filter(!is.na(MEDLYN_g1)) %>%
  select(SITE_NEON, LAYER_L, RAD_LW_ABS, RAD_SW_ABS) %>%
  mutate(LAYER_L = factor(LAYER_L)) %>%
  pivot_longer(c(RAD_SW_ABS, RAD_LW_ABS)) %>%
  ggplot(aes(x=factor(LAYER_L), y=value)) +
  geom_boxplot(aes(fill=name)) +
  facet_wrap(~ SITE_NEON, scales="free") +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  scale_fill_hue(labels=c("Longwave radiation", "Shortwave radiation")) +
  theme_bw() +
  labs(y=expression("Absorbed radiation"~bgroup("(", frac(W, m^2), ")")),
       x=expression("Cumulative LAI"~bgroup("(", frac(m^2, m^2), ")")),
       fill=""))

# Layer by layer plot of GPP
(plot_layer_gpp <- flux_ps_model %>%
  filter(!is.na(MEDLYN_g1)) %>%
  select(SITE_NEON, LAYER_L, PS_LAYER_GPP) %>%
  mutate(LAYER_L = factor(LAYER_L)) %>%
  ggplot(aes(x=factor(LAYER_L), y=PS_LAYER_GPP)) +
  geom_boxplot(fill="lightgreen") +
  facet_wrap(~ SITE_NEON, scales="free") +
  coord_flip() +
  scale_x_discrete(limits=rev)  +
  theme_bw() +
  labs(y=expression("Layer GPP"~bgroup("(", frac(mol~CO[2], m^2~s), ")")),
       x=expression("Cumulative LAI"~bgroup("(", frac(m^2, m^2), ")")),
       fill=""))

(plot_layer_gs <- flux_gs %>%
    filter(!is.na(PS_LAYER_GS)) %>%
    select(SITE_NEON, LAYER_L, PS_LAYER_GS) %>%
    mutate(LAYER_L = factor(LAYER_L)) %>%
    ggplot(aes(x=factor(LAYER_L), y=PS_LAYER_GS)) +
    geom_boxplot(fill="lightgreen") +
    facet_wrap(~ SITE_NEON, scales="free") +
    coord_flip() +
    scale_x_discrete(limits=rev)  +
    theme_bw() +
    labs(y=expression("Layer"~g[s]~bgroup("(", frac(mol, m^2~s), ")")),
         x=expression("Cumulative LAI"~bgroup("(", frac(m^2, m^2), ")")),
         fill=""))
