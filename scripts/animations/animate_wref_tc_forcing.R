library(tidyverse)
library(gganimate)
library(ggdark)
library(lidR)
source("scripts/energy_balance/leaf_eb_still.R")

theme_set(dark_theme_bw())

wref_tc <- read_csv("data_out/cross_site_gs_tc_result.csv") %>%
  filter(SITE_NEON == "WREF") %>%
  mutate(TIMESTAMP_LOCAL = TIMESTAMP - hours(8))

# Make sure timestamp looks right
ggplot(wref_tc, aes(x=hour(TIMESTAMP_LOCAL))) +
  geom_boxplot(aes(y=SW_IN_TOC, group=hour(TIMESTAMP_LOCAL)))

wref_a <- c(25, 8, 0.4, 20, 0.0725)

lidar_zmin <- 1
lidar_zmax <- 58

wref_lidar <- read_csv("data_out/cross_site_lidar_constants.csv") %>%
  filter(site == "WREF")

lad_profile <- readLAS("data_working/neon_neartower_lidar/WREF.laz") %>%
  filter_poi(Z < lidar_zmax, Z > lidar_zmin) %>%
  .$Z %>% LAD(k = wref_lidar$k_lai) %>%
  map_df(rev) %>%
  mutate(z = z + 0.5,
         lai = cumsum(lad),
         iv_io = exp(-wref_lidar$k_rad * lai))

ggplot(lad_profile, aes(x=iv_io, y=z)) +
  geom_line()

z_to_model <- lad_profile$z

# Interpolate weather conditions over the entire day
wref_one_day_weather <- wref_tc %>%
  filter(as.Date(TIMESTAMP_LOCAL) == ymd("2020-07-01")) %>%
  select(TIMESTAMP_LOCAL, z, 
         PA, SW_IN_TOC, PPFD_IN_TOC, LW_IN, LW_OUT, SW_OUT,
         ta_interp, ws_interp, par_interp, mr_interp) %>%
  group_by(TIMESTAMP_LOCAL) %>%
  nest(data=-c(TIMESTAMP_LOCAL, PA, SW_IN_TOC, SW_OUT, PPFD_IN_TOC, LW_IN, LW_OUT))

wref_z_interpolated <- wref_one_day_weather %>%
  mutate(
    data_allz = map(
      data, function(d) {
        d_allz <- data.frame(z=z_to_model)
        
        d_allz$ta_interp <- approx(d$z, d$ta_interp, z_to_model, rule=2)$y
        d_allz$ws_interp <- approx(d$z, d$ws_interp, z_to_model, rule=2)$y
        d_allz$par_interp <- approx(d$z, d$par_interp, z_to_model, rule=2)$y
        d_allz$mr_interp <- approx(d$z, d$mr_interp, z_to_model, rule=2)$y
        
        return(d_allz)
      }
    )
  ) %>%
  select(-data) %>%
  unnest(data_allz) %>%
  ungroup() %>% 
  # attach lidar metrics, calculate downwelling shortwave
  left_join(lad_profile, by="z") %>%
  mutate(sw_in = SW_IN_TOC * iv_io,
         sw_out = (SW_OUT / SW_IN_TOC) * sw_in) %>%
  # calculate condition of the air
  mutate(
    e    = PA * (1 + (mr_interp/1000)^-1)^-1,
    esat = saturation_vapor_pressure(ta_interp + 273.15),
    rh   = e / esat,
    vpd  = esat - e
  ) %>%
  # run gs model
  mutate(
    gs_interp = jarvis_gs(ta_interp, vpd, sw_in, wref_a) * 10
  )

# Run temperature model
gs_notransp <- 0
wref_tc_model <- wref_z_interpolated %>%
  mutate(
    G = 0, # ignore soil heat flux
    t_canopy_model_still = pmap(
      list(ta_interp, PA, rh, ws_interp, gs_interp, sw_in, sw_out, LW_IN, G),
      function(Ta, Pa, RH, u, gs, sw_in, sw_out, lw_in, G) {
        leaf_temperature_isothermal(
          Ta+273.15, Pa, RH, u, gs, sw_in, 0, sw_out, lw_in, G, eb_constants
        ) %>% as.data.frame()
      }
    ),
    t_notransp_model_still = pmap(
      list(ta_interp, PA, rh, ws_interp, sw_in, sw_out, LW_IN, G),
      function(Ta, Pa, RH, u, sw_in, sw_out, lw_in, G) {
        leaf_temperature_isothermal(
          Ta+273.15, Pa, RH, u, gs_notransp, sw_in, 0, sw_out, lw_in, G, eb_constants
        ) %>% as.data.frame()
      }
    )
  ) %>%
  unnest(t_canopy_model_still, names_sep="_") %>%
  unnest(t_notransp_model_still, names_sep="_") %>%
  mutate(
    # Convert temperature outputs into deg C
    t_canopy_model_still_Tl = t_canopy_model_still_Tl - 273.15,
    t_notransp_model_still_Tl = t_notransp_model_still_Tl - 273.15,
    # Calculate radiation and latent heat forcings
    dT = t_canopy_model_still_Tl - ta_interp,
    Rn_forcing = t_notransp_model_still_Tl - ta_interp,
    LE_forcing = t_canopy_model_still_Tl - t_notransp_model_still_Tl
  )
      
anim <- ggplot(wref_tc_model, aes(x=z)) +
  geom_line(aes(y=Rn_forcing, color="Rn")) +
  geom_line(aes(y=LE_forcing, color="LE")) +
  geom_line(aes(y=dT, color="T_leaf - T_air")) +
  coord_flip() +
  geom_hline(yintercept=0, color='white', linetype="dashed") +
  scale_color_manual(
    values=c(
      "T_leaf - T_air" = "#00BA38",
      "Rn" = "#F0E442",
      "LE" = "#56B4E9"
    )
  ) +
  ylim(-4, 4) +
  labs(y="Temperature forcing (K)",
       x="Canopy height (m)",
       title="{strftime(frame_time, format='%I:%M %p')}") +
  transition_time(TIMESTAMP_LOCAL)

anim
