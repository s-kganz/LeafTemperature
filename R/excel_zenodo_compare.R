library(tidyverse)
library(ggridges)
library(readxl)
library(ggpmisc)

excel <- read_xlsx(
  "data_in/chris_still_wref_tir/WR_TIRtemp,Met,Flux_1hr_140101-151231_160425-with-metadata-column-headers.xlsx",
  sheet=2,
  skip=1,
  col_types="numeric"
) %>%
  select(Yr, DOY, Hr, `S_in(W m-2)`, `Ta_70(oC)`, contains("Tcan_S")) %>%
  mutate(time_str = str_c(Yr, DOY, Hr, sep=" "),
         TIMESTAMP = strptime(time_str, format="%Y %j %H", tz="UTC"),
         TIMESTAMP = as_datetime(TIMESTAMP)) %>%
  rename(TA_excel = `Ta_70(oC)`)

zenodo <- read_csv("data_in/chris_still_wref_tir/still_WR_flux_data.csv") %>%
  select(TIMESTAMP_START, SW_IN, TA, contains("Tcan_S")) %>%
  mutate(
    TIMESTAMP_START = as.character(TIMESTAMP_START),
    TIMESTAMP_START = parse_datetime(TIMESTAMP_START, format="%Y%m%d%H%M")
  ) %>%
  rename(TA_zenodo = TA)

# Verify timestamps are in local time - peak shortwave flux at noon.
ggplot(excel, aes(x=hour(TIMESTAMP), y=`S_in(W m-2)`)) +
  geom_boxplot(aes(group=hour(TIMESTAMP)))

ggplot(zenodo, aes(x=hour(TIMESTAMP_START), y=SW_IN)) +
  geom_boxplot(aes(group=hour(TIMESTAMP_START)))

# Convert to long format and join by timestamp
excel_long <- excel %>%
  select(TIMESTAMP, TA_excel, contains("Tcan_S")) %>%
  pivot_longer(contains("Tcan_S"), names_to="roi", values_to="tcan_excel") %>%
  drop_na()

zenodo_long <- zenodo %>%
  select(TIMESTAMP_START, TA_zenodo, contains("Tcan_S")) %>%
  rename(TIMESTAMP = TIMESTAMP_START) %>%
  pivot_longer(contains("Tcan_S"), names_to="roi", values_to="tcan_zenodo") %>%
  drop_na()

both_long <- excel_long %>%
  inner_join(zenodo_long, 
             by=c("TIMESTAMP", "roi"))

# 1:1 plots
ggplot(both_long, aes(x=tcan_excel, y=tcan_zenodo)) +
  geom_point(alpha=0.5) +
  coord_equal() +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  facet_wrap(~ roi, nrow = 3) +
  theme_bw()

# Histogram of differences
ggplot(both_long, aes(x=tcan_excel - tcan_zenodo)) +
  geom_histogram() +
  geom_vline(xintercept=0, color="red", linetype="dashed") +
  facet_wrap(~ roi) +
  theme_bw() +
  labs(x="T difference (Excel - Zenodo)")

# TL - TA for both methods
both_long %>%
  mutate(
    dT_excel = tcan_excel - TA_excel,
    dT_zenodo = tcan_zenodo - TA_zenodo
  ) %>%
  select(TIMESTAMP, roi, dT_excel, dT_zenodo) %>%
  pivot_longer(contains("dT")) %>%
  ggplot() +
  geom_density_ridges(aes(x=value, y=name)) +
  scale_y_discrete(labels=c("Excel", "Zenodo")) +
  #facet_wrap(~ roi) +
  theme_bw() +
  labs(x="Tcan - TA (K)", y="")

# 1:1 of average Tcan across all ROIs
both_long %>%
  group_by(TIMESTAMP) %>%
  summarize(
    excel_mean_tcan = mean(tcan_excel),
    zenodo_mean_tcan = mean(tcan_zenodo)
  ) %>%
  ggplot(aes(x=excel_mean_tcan, y=zenodo_mean_tcan)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  theme_bw() +
  labs(x="Excel mean Tcan (C)",
       y="Zenodo mean Tcan (C)")

# Remake Fig. 2B
zenodo %>%
  mutate(
    TOC_Tcan = rowMeans(zenodo %>% select(matches("Tcan_S[0-4]")))
  ) %>%
  filter(month(TIMESTAMP_START) %in% 5:9,
         SW_IN > 25) %>%
  ggplot(aes(x=TA_zenodo, y=TOC_Tcan)) +
  geom_point(aes(color=hour(TIMESTAMP_START))) +
  stat_ma_line(linetype="dashed", color="red") +
  geom_abline(slope=1, intercept=0, color="red") +
  scale_color_viridis_c() +
  theme_bw() +
  labs(x="Air temp (C)",
       y="Tcan (C)",
       color="Hour of day")
  
             