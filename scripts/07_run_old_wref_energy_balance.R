# This script was written for use with the flux and TIR data from the Still 
# et al. PNAS paper. The source data for this is available at
# https://zenodo.org/records/6862565. The script is here in case anyone else
# wants to play with these data, but the analysis is not included in our paper.

library(tidyverse)
library(ggpmisc)
library(suncalc)
library(photosynthesis)
library(lmodel2)
source("scripts/energy_balance/leaf_eb_still.R")
source("scripts/energy_balance/leaf_eb_numeric.R")
source("scripts/tower_util.R")

# Read source data and constants ----

## Top-of-canopy meteorology ----
wref_tower <- read_csv(
  "data_out/old_wref_clean.csv",
  col_types=cols(TIMESTAMP=col_datetime(), .default=col_double()))

wref_toc <- wref_tower %>%
  select(TIMESTAMP, PA, P, G, SW_IN, SW_DIF, SW_OUT, LW_IN, LW_OUT, H, LE, 
         NEE, USTAR, TA_70, RH_70, WS_70, CO2, contains("Tcan")) %>%
  # Growing season only
  filter(month(TIMESTAMP) %in% 5:9)

# SW_DIF lookup table ----
# SW_DIF is completely NA where we have usable flux measurements. So instead we
# estimate the proportion of SW_DIF from a lookup table created with the NEON
# dataset.
neon_sw_dif <- read_csv("data_working/neon_flux/US-xWR.csv") %>%
  select(TIMESTAMP, matches("SW_IN_\\d_\\d_\\d"), SW_DIF)

sw_dif_lut <- neon_sw_dif %>%
  mutate(
    SW_IN = rowMeans(neon_sw_dif %>% select(matches("SW_IN_\\d_\\d_\\d")), 
                     na.rm=TRUE),
    PROP_DIF = pmax(1e-4, pmin(1, SW_DIF / SW_IN)),
    month = month(TIMESTAMP),
    hour = hour(TIMESTAMP)
  ) %>%
  filter(!is.na(PROP_DIF), month(TIMESTAMP) %in% 5:9) %>%
  group_by(month, hour) %>%
  summarize(
    # geometric mean since this is a ratio
    PROP_DIF = PROP_DIF %>% log() %>% mean() %>% exp()
  )

## Partitioned fluxes ----
wref_flux <- read_csv("data_out/old_wref_flux_partition.csv") %>%
  select(DateTime, NEE, Reco_uStar, GPP_uStar_f) %>%
  rename(TIMESTAMP=DateTime,
         TOC_NEE=NEE,
         TOC_Reco=Reco_uStar,
         TOC_GPP=GPP_uStar_f)

## Interpolated meteorology ----
wref_interp <- read_csv("data_out/old_wref_interpolated_meteorology.csv") %>%
  rename_with(
    .fn = function(x) str_c("LAYER_", x),
    .cols = -c(TIMESTAMP, layer_l)
  ) %>%
  rename(LAYER_L=layer_l)

## Constants & metadata ----
wref_meta <- read_csv("data_working/neon_site_metadata.csv") %>%
  filter(site_neon == "WREF")

aq_const <- read_csv("data_out/cross_site_aq_constants.csv") %>%
  filter(site == "WREF")

medlyn_roughness <- readRDS("data_out/old_wref_medlyn_roughness_params.rds")

lidar <- read_csv("data_out/neon_lidar_constants.csv") %>%
  filter(site == "WREF")

# Attach constants to flux table ----
wref_flux$lai           <- wref_meta$lai
# Calculate layer positions from LAI
wref_flux$LAYER_L       <- map(wref_flux$lai, function(x) seq(0, floor(x)))
wref_flux$TOWER_LAT     <- wref_meta$tower_lat
wref_flux$TOWER_LON     <- wref_meta$tower_lon
wref_flux$AQ_k_sat      <- aq_const$k_sat
wref_flux$AQ_phi_J      <- aq_const$phi_J
wref_flux$AQ_theta_J    <- aq_const$theta_J
wref_flux$MEDLYN_g0     <- medlyn_roughness$g0
wref_flux$MEDLYN_g1     <- medlyn_roughness$g1
# n.b. roughness parameters aren't needed in the all-NEON run because we have wind
# sensors in the canopy. Here we don't have enough, so we need to estimate
# wind speed based on the roughness pars.
wref_flux$ROUGHNESS_d   <- medlyn_roughness$d
wref_flux$ROUGHNESS_z0m <- medlyn_roughness$z0m
wref_flux$LIDAR_kb_real <- lidar$kb_real
wref_flux$LIDAR_kd_real <- lidar$kd_real

# Attach SW_DIF to toc table ----
wref_toc_sw_dif <- wref_toc %>%
  mutate(
    month = month(TIMESTAMP),
    hour = hour(TIMESTAMP)
  ) %>%
  left_join(sw_dif_lut, by=c("month", "hour")) %>%
  select(-month, -hour)

# Attach top-of-canopy rad to flux table ----
wref_toc_rad <- wref_toc_sw_dif %>%
  select(TIMESTAMP, SW_IN, PROP_DIF, SW_OUT, LW_IN, LW_OUT, H, LE, USTAR, PA, TA_70,
         WS_70, CO2) %>%
  rename_with(
    .fn = function(x) str_c("TOC_", x),
    .cols = -c(TIMESTAMP)
  ) %>%
  rename(TOC_TA=TOC_TA_70,
         TOC_WS=TOC_WS_70)

wref_flux_rad <- inner_join(
  wref_flux, wref_toc_rad, by="TIMESTAMP"
) %>%
  mutate(TOC_SW_DIF = TOC_PROP_DIF * TOC_SW_IN)

# Calculate sunlight position ----
wref_solar <- wref_flux_rad %>% bind_cols(
  getSunlightPosition(
    data = data.frame(
      date=with_tz(.$TIMESTAMP, tzone="Etc/GMT-8"),
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

# Expand layer L and attach within-canopy meteorology ----
wref_layer_l <- wref_solar %>%
  unnest(LAYER_L) %>%
  relocate(LAYER_L, .after=last_col()) %>%
  inner_join(
    wref_interp,
    by=c("TIMESTAMP"="TIMESTAMP", "LAYER_L")
  )

# Calculate wind speed based on top/bottom of canopy and roughness parameters. 
# The equation in bigleaf is only valid for the uppermost part of the canopy.
# Below z = z0m + d, just use the lower-canopy value.
z_critical <- medlyn_roughness$z0m + medlyn_roughness$d + 3 # buffer

wref_wind <- wref_layer_l %>%
  # Temporarily join the lower-canopy wind speed measurement
  left_join(
    wref_tower %>% select(TIMESTAMP, WS_02),
    by="TIMESTAMP"
  ) %>%
  mutate(
    LAYER_WS = ifelse(
      LAYER_z < z_critical,
      WS_02,
      bigleaf::wind.profile(
          data=NULL, z=LAYER_z, Tair=.$TOC_TA, pressure=.$TOC_PA,
          ustar=.$TOC_USTAR, H=.$TOC_H, wind=.$TOC_WS,
          zr=wref_meta$tower_height, zh=wref_meta$canopy_height,
          d=medlyn_roughness$d, z0m=medlyn_roughness$z0m
        )
    ),
    # Ensure nonzero windspeed
    LAYER_WS = ifelse(LAYER_WS == 0, WS_02, LAYER_WS)
  ) %>%
  select(-WS_02)

# Drop nighttime observations, those outside the growing season, and NAs
wref_filter <- wref_wind %>%
  filter(TOC_SW_IN > 100,
         month(TIMESTAMP) %in% 5:9) %>%
  drop_na()

# Run radiation submodel ----
flux_rad_model <- wref_filter %>%
  mutate(
    RAD_SKY_EMIS = TOC_LW_IN / (eb_constants_$sb * (TOC_TA+273.15)^4),
    RAD_SKY_EMIS = pmax(0, pmin(1, RAD_SKY_EMIS)),
    RAD_SW_ABS = absorbed_shortwave_radiation(
      SW_beam=pmax(0, TOC_SW_IN - TOC_SW_DIF),
      SW_dif=TOC_SW_DIF,
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
      Ta=(TOC_TA+273.15),
      l=LAYER_L,
      L=lai,
      kb=-LIDAR_kb_black,
      ea=RAD_SKY_EMIS,
      ef=0.95,
      es=0.9,
      kd=-LIDAR_kd_black
    )
  )

# Run photosynthesis submodel ----
flux_ps_weight <- flux_rad_model %>%
  # Calculate layer contribution to GPP
  mutate(
    PS_LAYER_GPP_WEIGHT = marshall_biscoe_1980(
      RAD_SW_ABS, AQ_k_sat, AQ_phi_J, AQ_theta_J
    )
  ) %>%
  # Sum the contributions within the whole canopy
  group_by(TIMESTAMP) %>%
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
    LAYER_VPD = vapor_pressure_deficit(LAYER_TA+273.15, LAYER_RH),
    LAYER_VP = LAYER_SAT_VP - LAYER_VPD,
  )

flux_ps_model <- flux_ps_weight %>%
  # Calculate stomatal conductance in the layer
  mutate(PS_LAYER_GS = medlyn_gs(PS_LAYER_GPP, LAYER_VPD, TOC_CO2, MEDLYN_g1, MEDLYN_g0))

# Run the model ----
flux_eb_result <- flux_ps_model %>%
  mutate(
    EB_MODEL = pmap(
      list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, TOC_PA, LAYER_RH,
           RAD_SW_ABS, RAD_LW_ABS),
      energy_balance_driver,
      bounds=20,
      run_no_transp=TRUE,
      .progress="numeric EB model"
    ),
    # EB_MODEL = pmap(
    #   list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, TOC_PA, LAYER_RH,
    #        RAD_SW_ABS, RAD_LW_ABS),
    #   leaf_temperature_isothermal,
    #   aslist=TRUE,
    #   .progress="isothermal EB model"
    # )
  ) %>%
  # unnest_wider(EB_MODEL_ISOTHERMAL, names_sep="_") %>%
  unnest_wider(EB_MODEL, names_sep="_") %>%
  mutate(EB_MODEL_dT = EB_MODEL_dT_Rn - EB_MODEL_dT_LE)

# Save results
write_if_not_exist(flux_eb_result, "data_out/model_runs/old_wref_eb.csv")

# Validate against the thermal camera measurements ----
tir_heights <- data.frame(
  field=str_c("Tcan_S", 0:9),
  zmin=c(rep(50, 5), rep(15, 3), rep(2, 2)),
  zmax=c(rep(80, 5), rep(40, 3), rep(10, 2)),
  group=c(rep("Upper canopy", 5), rep("Middle canopy", 3), rep("Lower canopy", 2))
)

wref_tir_Tl <- wref_tower %>%
  select(TIMESTAMP, contains("Tcan"), -Tcan_avg) %>%
  pivot_longer(-c(TIMESTAMP)) %>%
  left_join(
    tir_heights,
    by=c("name"="field")
  ) %>%
  filter(!is.na(value)) %>%
  group_by(TIMESTAMP, zmin, zmax, group) %>%
  summarize(TIR_Tl = mean(value))

wref_model_Tl <- flux_eb_result %>%
  mutate(EB_MODEL_Tl = EB_MODEL_Tl - 273.15) %>%
  select(TIMESTAMP, LAYER_z, LAYER_TA, EB_MODEL_Tl) %>%
  left_join(
    wref_tir_Tl,
    join_by(TIMESTAMP == TIMESTAMP, LAYER_z >= zmin, LAYER_z <= zmax)
  )
  group_by(TIMESTAMP, group) %>%
  summarize(
    across(c(EB_MODEL_Tl, LAYER_TA), mean)
  )

# Join in the Ameriflux air temperature data because it is better
amf_wrc <- amerifluxr::amf_read_base(
  file="data_in/AMF_US-Wrc_BASE-BADM_8-1/AMF_US-Wrc_BASE_HH_8-1.csv",
  unzip=FALSE,
  parse_timestamp=TRUE
) %>%
  mutate(TIMESTAMP = TIMESTAMP - minutes(15) - hours(8))

amf_wrc_filter <- amf_wrc %>%
  filter(TIMESTAMP > min(wref_tower$TIMESTAMP),
         TIMESTAMP < max(wref_tower$TIMESTAMP))

wref_model_Tl_amf_ta <- wref_model_Tl %>%
  inner_join(
    amf_wrc_filter %>% 
      select(TIMESTAMP, TA) %>%
      rename(TA_AMF = TA),
    by="TIMESTAMP"
  )

wref_model_Tl_amf_ta %>%
  select(TA_AMF, group, TIR_Tl, EB_MODEL_Tl) %>%
  pivot_longer(-c(TA_AMF, group)) %>%
  mutate(name = recode(
    name,
    "TIR_Tl"="Thermal IR",
    "EB_MODEL_Tl"="EB Model"
  )) %>%
  filter(name == "Thermal IR") %>%
  ggplot(aes(x=TA_AMF, y=value, group=name)) + 
  geom_point(alpha=0.4, aes(color=name)) +
  stat_ma_line(se=FALSE, color="black", show.legend=FALSE) +
  stat_ma_eq(use_label(c("grp", "eq"))) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  facet_wrap(~ factor(group, levels=c("Upper canopy", "Middle canopy", "Lower canopy")), 
             nrow=1, strip.position="top") +
  coord_equal() +
  labs(x=expression("Air temperature"~"("*degree*"C)"),
       y=expression("Leaf temperature"~"("*degree*"C)"),
       color="") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0),
        strip.background = element_rect(fill=NA),
        strip.text.x = element_text(size=12),
        legend.text = element_text(size=12)) +
  ggtitle("Ameriflux TA")
  

amf_ta_join <- amf_wrc_filter %>%
  select(TIMESTAMP, TA) %>%
  rename(TA_amf = TA) %>%
  inner_join(
    wref_tower %>%
      select(TIMESTAMP, TA_70) %>%
      rename(TA_excel=TA_70),
    by="TIMESTAMP"
  )

ggplot(amf_ta_join, aes(x=TA_amf, y=TA_excel)) +
  geom_hex() +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  coord_fixed() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x=expression("Ameriflux TA"~"("*degree*"C)"),
       y=expression("Excel sheet TA"~"("*degree*"C)"))

       