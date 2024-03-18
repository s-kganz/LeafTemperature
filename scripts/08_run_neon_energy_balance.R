# This script runs the final energy balance model.

library(tidyverse)
library(ggpointdensity)
library(ggridges)
library(suncalc)
library(photosynthesis)
library(lmodel2)
source("scripts/energy_balance/leaf_eb_still.R")
source("scripts/energy_balance/leaf_eb_numeric.R")
source("scripts/tower_util.R")
# Read source data and constants ----

## Whole-canopy fluxes ----
flux <- read_csv("data_out/cross_site_flux_partition_qc.csv") %>%
  select(site, timeBgn, GPP_uStar_f, data.fluxH2o.nsae.flux) %>%
  rename(SITE_NEON=site,
         TIMESTAMP=timeBgn,
         TOC_GPP=GPP_uStar_f,
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
    LAYER_L = map(LAI, function(x) seq(0, floor(x)))
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
flux_ps_weight <- flux_rad_model %>%
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
  )
  
flux_ps_model <- flux_ps_weight %>%
  # Drop observations with missing data for getting stomatal conductance
  drop_na(c(MEDLYN_g1, LAYER_CO2, LAYER_H2O)) %>%
  # Calculate stomatal conductance in the layer
  mutate(PS_LAYER_GS = medlyn_gs(PS_LAYER_GPP, LAYER_VPD, LAYER_CO2, MEDLYN_g1, MEDLYN_g0))

# Run energy balance models ----
flux_eb_result <- flux_ps_model %>%
  mutate(
    EB_MODEL = pmap(
      list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, RAD_TOC_PA, LAYER_RH,
           RAD_SW_ABS, RAD_LW_ABS),
      energy_balance_driver,
      bounds=20,
      run_no_transp=TRUE,
      .progress="numeric EB model"
    ),
    # EB_MODEL_ISOTHERMAL = pmap(
    #   list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, RAD_TOC_PA, LAYER_RH,
    #        RAD_SW_ABS, RAD_LW_ABS),
    #   leaf_temperature_isothermal,
    #   aslist=TRUE,
    #   .progress="isothermal EB model"
    # ),
  ) %>%
  # unnest_wider(EB_MODEL_ISOTHERMAL, names_sep="_") %>%
  unnest_wider(EB_MODEL, names_sep="_")
  

# Second run with gbH way smaller by increasing characteristic leaf dimension
# by an order of magnitude.
flux_eb_result_bigleaf <- flux_ps_model %>%
  mutate(
    EB_MODEL = pmap(
      list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, RAD_TOC_PA, LAYER_RH,
           RAD_SW_ABS, RAD_LW_ABS),
      energy_balance_driver,
      d=0.1,
      bounds=20,
      run_no_transp=TRUE,
      .progress="low gbH EB model"
    )
  ) %>%
  unnest_wider(EB_MODEL, names_sep="_") %>%
  mutate(EB_MODEL_dT = EB_MODEL_Tl - 273.15 - LAYER_TA)
  

# Third run on WREF only, and varying g1 from 2-10 for sensitivity analysis.
flux_wref_g1_sensitivity <- flux_ps_weight %>%
  filter(SITE_NEON == "WREF") %>%
  drop_na(c(LAYER_H2O, LAYER_CO2)) %>%
  mutate(MEDLYN_g0 = 0, MEDLYN_g1 = map(SITE_NEON, function(x) seq(2:10))) %>%
  unnest(MEDLYN_g1) %>%
  mutate(PS_LAYER_GS = medlyn_gs(PS_LAYER_GPP, LAYER_VPD, LAYER_CO2, MEDLYN_g1, MEDLYN_g0))

flux_eb_result_wref_g1_sensitivity <- flux_wref_g1_sensitivity %>%
  mutate(
    EB_MODEL = pmap(
      list(LAYER_TA+273.15, LAYER_WS, PS_LAYER_GS, RAD_TOC_PA, LAYER_RH,
           RAD_SW_ABS, RAD_LW_ABS),
      energy_balance_driver,
      bounds=20,
      run_no_transp=TRUE,
      .progress="WREF g1 sensitivty analysis"
    )
  ) %>%
  unnest_wider(EB_MODEL, names_sep="_") %>%
  mutate(EB_MODEL_dT = EB_MODEL_Tl - 273.15 - LAYER_TA)
  
# Save results ----
write_if_not_exist(flux_eb_result, "data_out/model_runs/cross_site_eb.csv")
write_if_not_exist(flux_eb_result_bigleaf, "data_out/model_runs/cross_site_eb_bigleaf.csv")
write_if_not_exist(flux_eb_result_wref_g1_sensitivity, "data_out/model_runs/wref_eb_g1_sensitivity.csv")

# Plots! ----
# Tl by site
(plot_tl_hist <- flux_eb_result_bigleaf %>%
   #filter(LAYER_L/LAI > 0.5) %>%
   ggplot(aes(x=LAYER_TA, y=EB_MODEL_Tl-273.15-LAYER_TA)) +
   geom_point(aes(color=LAYER_L/LAI), alpha=0.5) +
   #geom_abline(slope=1, intercept=0, color="red") +
   facet_wrap(~ SITE_NEON) +
   scale_color_viridis_c() +
   theme_bw()
)

# Tl vs. Ta for all layers
(plot_tl_ta <- flux_eb_result %>%
  mutate(EB_MODEL_dT = EB_MODEL_Tl - 273.15 - LAYER_TA) %>%
  ggplot(aes(x=LAYER_TA, y=EB_MODEL_Tl-273.15)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=0) +
  facet_grid(SITE_NEON ~ LAYER_L) +
  theme_bw() +
  labs(x=expression("Air temperature ("~degree~"C)"),
       y="Leaf-air T difference ("~degree~"K)")
)

# Four quadrants figure
flux_eb_models <- flux_eb_result %>%
  mutate(EB_MODEL_dT = EB_MODEL_Tl - 273.15 - LAYER_TA) %>%
  select(SITE_NEON, LAI, LAYER_L, LAYER_TA, contains("EB_MODEL")) %>%
  group_by(SITE_NEON, LAYER_L, LAI) %>%
  nest(data=-c(SITE_NEON, LAYER_L, LAI)) %>%
  mutate(
    # Run model II regressions to get Tl/Ta slope
    slope_lm = map(
      data,
      function(data) {
        lmodel2(EB_MODEL_Tl ~ LAYER_TA, data=data)
      }
    ),
    slope_p025 = map_dbl(slope_lm, function(m) m$confidence.intervals[3, 4]),
    slope_p50  = map_dbl(slope_lm, function(m) m$regression.results[3, 3]),
    slope_p975 = map_dbl(slope_lm, function(m) m$confidence.intervals[3, 5]),
    int_p50    = map_dbl(slope_lm, function(m) m$regression.results[3, 2] - 273.15),
    # Variance partitioning
    var_lm = map(
      data,
      function(data) {
        anova(lm(EB_MODEL_dT ~ EB_MODEL_dT_Rn + EB_MODEL_dT_LE, data=data))
      }
    ),
    tot_var = map_dbl(var_lm, function(m) sum(m$`Sum Sq`)),
    Rn_var = map_dbl(var_lm, function(m) m$`Sum Sq`[1]),
    LE_var = map_dbl(var_lm, function(m) sum(m$`Sum Sq`[2])),
    prop_var_Rn = Rn_var / tot_var,
    prop_var_LE = LE_var / tot_var
  ) %>%
  select(SITE_NEON, LAYER_L, LAI, prop_var_Rn, slope_p50, slope_p025, slope_p975, int_p50) %>%
  arrange(SITE_NEON, LAYER_L) %>%
  group_by(SITE_NEON) %>%
  mutate(
    xend = lead(prop_var_Rn),
    yend = lead(slope_p50)
  ) 

(plot_slope_rn <- flux_eb_models %>%
  ggplot(aes(x=prop_var_Rn, y=slope_p50)) +
  #geom_segment(aes(group=SITE_NEON, xend=xend, yend=yend)) +
  geom_errorbar(aes(ymin=slope_p025, ymax=slope_p975),
               width=0.001) +
  geom_point(aes(color=LAYER_L/LAI), size=4) +
  scale_color_viridis_c() +
  facet_wrap(~ SITE_NEON, scales="free") +
  labs(y=expression(T[leaf]~"vs."~T[air]~"slope (model II regression)"),
       x=expression("Proportion of variance explained by"~R[n]),
       color="Proportion of total LAI") +
  theme_bw())

site_set_1 <- c("DEJU", "ABBY", "RMNP", "WREF")

(plot_regression_result <- flux_eb_models %>%
    filter(!(SITE_NEON %in% site_set_1)) %>%
    mutate(LAYER_L = factor(LAYER_L)) %>%
    ggplot(aes(x=slope_p50, y=LAYER_L)) +
    geom_vline(xintercept=1, color="grey50", linetype="dashed") +
    geom_errorbar(aes(xmin=slope_p025, xmax=slope_p975)) +
    geom_point(aes(fill=int_p50), color="black", size=4, pch=21) +
    geom_label(aes(label=SITE_NEON), x=-Inf, y=Inf,
               hjust="left", vjust="top", size=5) +
    scale_fill_distiller(
      type="div", palette="Spectral",
      limits=c(-1.5, 1.5)
    ) +
    scale_y_discrete(limits=rev) +
    scale_x_continuous(limits=c(0.97, 1.15)) +
    facet_grid(SITE_NEON ~ .,
               scales="free_y", space="free") +
    theme_bw() +
    theme(strip.placement="outside",
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          axis.title = element_text(size=18),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14)
          ) +
    guides(fill = guide_colorbar(
      ticks.colour="black",
      frame.colour="black"
    )) +
    labs(x=expression("T"[leaf]~"vs"~"T"[air]~"regression slope ("~degree*C*"/"*degree*C*")"),
         y=expression("Cumulative LAI (m"^2~"m"^-2*")"),
         fill=expression("Model II\nregression intercept"~"("*degree*C*")"))
)

# How well does layer-by-layer LE match canopy LE?
(plot_le_one_to_one <- flux_eb_result %>%
  select(SITE_NEON, TIMESTAMP, TOC_LE, EB_MODEL_LE) %>%
  filter(TOC_LE > 0) %>%
  group_by(SITE_NEON, TIMESTAMP) %>%
  summarize(TOC_LE=first(TOC_LE),
            EB_MODEL_LE_SUM = sum(EB_MODEL_LE)) %>%
  ggplot(aes(x=TOC_LE, y=EB_MODEL_LE_SUM)) +
  ggpointdensity::geom_pointdensity() +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  coord_equal() +
  scale_color_viridis_c() +
  facet_wrap(~ SITE_NEON) +
  theme_bw() +
  labs(x=expression("Tower-based "~lambda~"E"~bgroup("(",frac(W, m^2),")")),
       y=expression("Model-based "~lambda~"E"~bgroup("(",frac(W, m^2),")")),
       color="Point density")
)

# Layer-by-layer plot of the two fluxes
(plot_rad_fluxes <- flux_eb_result %>%
  filter(!is.na(MEDLYN_g1), SITE_NEON == "WREF") %>%
  select(SITE_NEON, LAYER_L, RAD_LW_ABS, RAD_SW_ABS) %>%
  mutate(LAYER_L = factor(LAYER_L)) %>%
  pivot_longer(c(RAD_SW_ABS, RAD_LW_ABS)) %>%
  ggplot(aes(x=factor(LAYER_L), y=value)) +
  geom_boxplot(aes(fill=name)) +
  facet_wrap(~ SITE_NEON, scales="free_y") +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  scale_fill_hue(labels=c("Longwave radiation", "Shortwave radiation")) +
  theme_bw() +
  labs(y=expression("Absorbed radiation (W m"^-2*")"),
       x=expression("Cumulative LAI (m"^2~"m"^-2*")"),
       fill="")) +
  theme(legend.text=element_text(size=18),
        axis.title=element_text(size=18))

# Layer by layer plot of GPP
(plot_layer_gpp <- flux_ps_model %>%
  filter(!is.na(MEDLYN_g1), SITE_NEON=="WREF") %>%
  select(SITE_NEON, LAYER_L, PS_LAYER_GPP) %>%
  mutate(LAYER_L = factor(LAYER_L)) %>%
  ggplot(aes(x=factor(LAYER_L), y=PS_LAYER_GPP)) +
  geom_boxplot(fill="lightgreen") +
  facet_wrap(~ SITE_NEON, scales="free") +
  coord_flip() +
  scale_x_discrete(limits=rev)  +
  theme_bw() +
  labs(y=expression("Layer GPP ("*mu*mol~CO[2]~m^-2~s^-1*")"),
       x=expression("Cumulative LAI (m"^2~"m"^-2*")"),
       fill="") +
  theme(axis.title = element_text(size=18)))

(plot_layer_gs <- flux_ps_model %>%
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

# Changes in Tl/Ta slope by varying g1 in WREF
(plot_wref_g1_tl_ta <- flux_eb_result_wref_g1_sensitivity %>%
    mutate(EB_MODEL_dT = EB_MODEL_Tl - 273.15 - LAYER_TA) %>%
    ggplot(aes(x=LAYER_TA, y=EB_MODEL_dT)) +
    geom_point(alpha=0.5) +
    geom_hline(yintercept=0) +
    facet_grid(MEDLYN_g1 ~ LAYER_L) +
    theme_bw() +
    labs(x=expression("Air temperature ("~degree~"C)"),
         y="Leaf-air T difference ("~degree~"K)")
)

# TOC gs by varying g1 in WREF
wref_g1 <- flux_eb_result %>% filter(SITE_NEON == "WREF") %>%
  pull(MEDLYN_g1) %>% head(1)
  
(plot_wref_g1_gs <- flux_eb_result_wref_g1_sensitivity %>%
    filter(LAYER_L == 0, PS_LAYER_GS < 0.5) %>%
    ggplot(aes(x=PS_LAYER_GS, y=MEDLYN_g1)) +
    geom_density_ridges(aes(group=MEDLYN_g1)) +
    geom_hline(yintercept=wref_g1, color="red", linetype="dashed") +
    scale_y_continuous(breaks=seq(1:9)) +
    theme_bw() +
    labs(x=expression("Layer"~g[s]~bgroup("(", frac(mol, m^2~s), ")")),
         y=expression("Medlyn"~g[1]~"parameter"))
)

# But, is total LE reasonable at high values of g1?
(plot_wref_g1_le_one_to_one <- flux_eb_result_wref_g1_sensitivity %>%
    select(SITE_NEON, TIMESTAMP, TOC_LE, EB_MODEL_LE, MEDLYN_g1) %>%
    filter(TOC_LE > 0) %>%
    mutate(MEDLYN_g1 = str_c("g1 = ", MEDLYN_g1)) %>%
    group_by(SITE_NEON, TIMESTAMP, MEDLYN_g1) %>%
    summarize(TOC_LE=first(TOC_LE),
              EB_MODEL_LE_SUM = sum(EB_MODEL_LE)) %>%
    ggplot(aes(x=TOC_LE, y=EB_MODEL_LE_SUM)) +
    ggpointdensity::geom_pointdensity() +
    geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
    coord_equal() +
    scale_color_viridis_c() +
    scale_x_continuous(breaks=seq(0, 450, by=150)) +
    facet_wrap(~ MEDLYN_g1) +
    theme_bw() +
    labs(x=expression("Tower-based "~lambda~"E"~bgroup("(",frac(W, m^2),")")),
         y=expression("Model-based "~lambda~"E"~bgroup("(",frac(W, m^2),")")),
         color="Point density")
)

# Tl vs. Ta slopes for different values of g1
wref_g1_ta_tl_slopes <- flux_eb_result_wref_g1_sensitivity %>%
  select(SITE_NEON, LAI, LAYER_L, MEDLYN_g1, LAYER_TA, contains("EB_MODEL")) %>%
  group_by(SITE_NEON, MEDLYN_g1, LAYER_L, LAI) %>%
  nest(data=-c(SITE_NEON, MEDLYN_g1, LAYER_L, LAI)) %>%
  mutate(
    # Run model II regressions to get Tl/Ta slope
    slope_lm = map(
      data,
      function(data) {
        lmodel2(EB_MODEL_Tl ~ LAYER_TA, data=data)
      }
    ),
    slope_p025 = map_dbl(slope_lm, function(m) m$confidence.intervals[3, 4]),
    slope_p50  = map_dbl(slope_lm, function(m) m$regression.results[3, 3]),
    slope_p975 = map_dbl(slope_lm, function(m) m$confidence.intervals[3, 5])) %>%
  ungroup()


# Compare relative flux magnitudes
flux_eb_result %>%
  filter(SITE_NEON == "WREF") %>%
  select(SITE_NEON, LAYER_L, EB_MODEL_Rn, EB_MODEL_H, EB_MODEL_LE) %>%
  pivot_longer(contains("EB_MODEL")) %>%
  ggplot(aes(x=value)) +
  geom_boxplot(aes(fill=name), position=position_dodge()) +
  facet_wrap(~ factor(LAYER_L), ncol=1, strip.position="left") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_blank()) +
  labs(x=expression("Heat flux"~"["~W/m^2~"]"),
       y="Cumulative LAI",
       fill="Flux source")
  
# TL vs. TA for just WREF
flux_eb_result %>%
  filter(SITE_NEON == "WREF") %>%
  mutate(canopy_pos = case_when(
    LAYER_L %in% 0:2 ~ "Upper canopy",
    LAYER_L %in% 3:6 ~ "Middle canopy",
    LAYER_L %in% 7:9 ~ "Lower canopy"
  )) %>%
  ggplot(aes(x=LAYER_TA, y=EB_MODEL_Tl-273.15)) +
  geom_abline(slope=1, intercept=0, color="red") +
  geom_point(alpha=0.5) +
  coord_equal() +
  facet_wrap(~ canopy_pos) +
  theme_bw() +
  theme(axis.title=element_text(size=18),
        strip.text=element_text(size=18)) +
  labs(x=expression("Air temperature ("*degree*"C)"),
       y=expression("Leaf temperature ("*degree*"C)"))

# Shaded GPP calculation
flux_eb_result %>%
  mutate(pct_sunlit = exp(LIDAR_kd_black * LAYER_L),
         GPP_sunlit = PS_LAYER_GPP * pct_sunlit,
         GPP_shade  = PS_LAYER_GPP * (1 - pct_sunlit)) %>%
  group_by(SITE_NEON) %>%
  summarize(tot_gpp = sum(PS_LAYER_GPP),
            tot_gpp_sunlit = sum(GPP_sunlit),
            tot_gpp_shade = sum(GPP_shade),
            prop_gpp_sunlit = tot_gpp_sunlit / tot_gpp,
            prop_gpp_shade = tot_gpp_shade / tot_gpp)
