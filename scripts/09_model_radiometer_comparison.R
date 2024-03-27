# Determine radiometer heights for each site
# Collect model layers for each site
# Figure out a way to match them up
# Make 1:1 plots of radiometer T vs. model leaf T. One per site in main text,
# then upper/lower canopy in supplemental?

library(tidyverse)
library(fuzzyjoin)
library(ggpointdensity)
source("scripts/tower_util.R")

# Read data ----
eb_result <- read_csv("data_out/model_runs/cross_site_eb.csv")
rad_tcan <- read_csv("data_out/cross_site_tcan.csv")

## Collect heights for model/radiometers ----
rad_heights <- amf_var_heights_ %>%
  filter(Site_ID %in% unique(eb_result$SITE_AMF),
         str_detect(Variable, "T_CANOPY_\\d_\\d_\\d"),
         Height > 2) %>%
  select(Site_ID, Variable, Height) %>%
  group_by(Site_ID) %>%
  mutate(Height_rank = rank(Height))

model_heights <- eb_result %>% select(SITE_AMF, LAYER_z) %>% unique()

rad_model_height_plot <- rad_heights %>%
  select(-Variable, -Height_rank) %>%
  rename(SITE_AMF=Site_ID,
         LAYER_z=Height) %>%
  mutate(instrument="Radiometer") %>%
  bind_rows(model_heights %>% mutate(instrument="Model")) %>%
  ggplot(aes(x=instrument, y=LAYER_z)) +
  geom_point() +
  facet_wrap(~ SITE_AMF, scales="free_y") +
  theme_bw()

# Pair up each model layer with the nearest radiometer ----
model_rad_join <- model_heights %>%
  # Cartesian product of radiometers/layers in each site
  full_join(rad_heights, by=c("SITE_AMF"="Site_ID"),
            relationship="many-to-many") %>%
  # Only keep the matches that minimize the difference in height
  # between the layer and its paired radiometer.
  mutate(absdist = abs(LAYER_z - Height)) %>%
  group_by(SITE_AMF, LAYER_z) %>%
  filter(absdist == min(absdist)) %>%
  ungroup() %>%
  select(-absdist)

rad_model_height_plot +
  geom_segment(data=model_rad_join, 
               mapping=aes(x="Model", xend="Radiometer", y=LAYER_z, yend=Height),
               color="red", linetype="dashed") +
  labs(x="", y="Height above ground (m)")

# Join radiometer Tcan with model output ----
## First join layer z to the tcan data
rad_tcan_layer_z <- rad_tcan %>%
  inner_join(model_rad_join %>% select(SITE_AMF, LAYER_z, Variable),
             by=c("SITE"="SITE_AMF", "TCAN_VAR"="Variable"))

eb_result_rad_tcan <- eb_result %>%
  inner_join(rad_tcan_layer_z, by=c("TIMESTAMP", "LAYER_z", "SITE_AMF"="SITE"))

# Summarize results for each pairing ----
eb_result_rad_tcan_summary <- eb_result_rad_tcan %>%
  group_by(TIMESTAMP, SITE_NEON, TCAN_VAR) %>%
  summarize(TCAN=mean(TCAN, na.rm=TRUE),
            LAYER_TA=mean(LAYER_TA, na.rm=TRUE),
            EB_MODEL_Tl=mean(EB_MODEL_Tl, na.rm=TRUE)-273.15)

write_if_not_exist(eb_result_rad_tcan_summary, "data_out/cross_site_model_rad_comparison.csv")
write_if_not_exist(model_rad_join, "data_out/model_rad_height_crosswalk.csv")

## Accuracy metrics ----
# RMSE, Bias, R2
accuracy_metrics <- eb_result_rad_tcan_summary %>%
  group_by(SITE_NEON) %>%
  summarize(
    rmse = sqrt(mean((TCAN - EB_MODEL_Tl)^2)),
    bias = mean(EB_MODEL_Tl - TCAN),
    r2   = cor(EB_MODEL_Tl, TCAN)^2
  )
# Plotting ----
ggplot(eb_result_rad_tcan_summary, aes(x=TCAN, y=EB_MODEL_Tl)) +
  geom_bin2d() +
  geom_abline(slope=1, intercept=0, color="black", linetype="dashed") +
  scale_fill_viridis_c(limits=c(0, 200), oob=scales::squish) +
  facet_wrap(~ SITE_NEON) +
  theme_bw() +
  coord_equal() +
  labs(x=expression("Radiometer temperature ("*degree*"C)"),
       y=expression("Model temperature ("*degree*"C)"),
       fill="Count")
