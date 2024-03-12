library(tidyverse)
library(vegan)
library(crossmatch)
library(foreach)
library(AICcPermanova)
library(lmodel2)
library(ggridges)
source("scripts/tower_util.R")

# Objective 1: validate EB model ----
# Two validation datasets: Chris's TIR data (2014-2015) and NEON radiometers
# installed after 2018.

## Join old WREB model output with TIR ----
old_wref_tower <- read_csv("data_out/old_wref_clean.csv")

tir_heights <- data.frame(
  field=str_c("Tcan_S", 0:9),
  zmin=c(rep(50, 5), rep(15, 3), rep(2, 2)),
  zmax=c(rep(80, 5), rep(40, 3), rep(10, 2)),
  group=c(rep("Upper canopy", 5), rep("Middle canopy", 3), rep("Lower canopy", 2))
)

wref_tir_Tl <- old_wref_tower %>%
  select(TIMESTAMP, contains("Tcan"), -Tcan_avg) %>%
  pivot_longer(-c(TIMESTAMP)) %>%
  left_join(
    tir_heights,
    by=c("name"="field")
  ) %>%
  filter(!is.na(value)) %>%
  group_by(TIMESTAMP, zmin, zmax, group) %>%
  summarize(TIR_Tl = mean(value))

old_wref_eb <- read_csv("data_out/model_runs/old_wref_eb.csv")

old_wref_model_tir <- old_wref_eb %>%
  mutate(EB_MODEL_Tl = EB_MODEL_Tl - 273.15) %>%
  select(TIMESTAMP, LAYER_z, LAYER_TA, EB_MODEL_Tl) %>%
  # Add TIR heights and summarize runs within each layer
  left_join(tir_heights, join_by(LAYER_z >= zmin, LAYER_z <= zmax)) %>%
  group_by(TIMESTAMP, group) %>%
  summarize(across(c(LAYER_TA, EB_MODEL_Tl), mean)) %>%
  inner_join(
    wref_tir_Tl,
    by=c("TIMESTAMP", "group")
  )

old_wref_mantel_table <- old_wref_model_tir %>%
  mutate(EB_MODEL_dT = EB_MODEL_Tl - LAYER_TA,
         TIR_dT = TIR_Tl - LAYER_TA) %>%
  select(TIMESTAMP, group, EB_MODEL_dT, TIR_dT) %>%
  pivot_longer(c(EB_MODEL_dT, TIR_dT)) %>%
  pivot_wider(names_from="group", values_from="value") %>%
  drop_na() %>%
  rename(temp_src = name)

## Join post-2018 model output with canopy radiometers
wref_tcan_height <- amf_var_heights_ %>%
  filter(Site_ID == "US-xWR") %>%
  filter(str_detect(Variable, "T_CANOPY")) %>%
  select(Variable, Height)

wref_rad <- read_csv("data_working/neon_flux/US-xWR.csv") %>%
  select(TIMESTAMP, contains("T_CANOPY")) %>%
  pivot_longer(-TIMESTAMP) %>%
  left_join(wref_tcan_height, by=c("name"="Variable")) %>%
  drop_na()

wref_eb <- read_csv("data_out/model_runs/cross_site_eb.csv") %>%
  filter(SITE_NEON == "WREF")

# Pair up model layers with corresponding radiometer and summarize
wref_eb_summary <- wref_eb %>%
  mutate(
    EB_MODEL_dT = EB_MODEL_Tl - LAYER_TA - 273.15,
    group = case_when(
      floor(LAYER_z) == 57 ~ 0,
      floor(LAYER_z) == 36 ~ 1,
      floor(LAYER_z) == 28 ~ 2,
      floor(LAYER_z) == 21 ~ 2,
      floor(LAYER_z) == 14 ~ 3,
      .default=NA
    )
  ) %>% filter(!is.na(group)) %>%
  group_by(TIMESTAMP, group) %>%
  summarize(EB_MODEL_dT = mean(EB_MODEL_dT),
            LAYER_TA = mean(LAYER_TA))

wref_rad_select <- wref_rad %>%
  mutate(
    group = case_when(
      floor(Height) == 53 ~ 0,
      floor(Height) == 35 ~ 1,
      floor(Height) == 23 ~ 2,
      floor(Height) == 12 ~ 3,
      .default=NA
    )
  ) %>% filter(!is.na(group)) %>%
  select(TIMESTAMP, value, group) %>%
  rename(RAD_Tl = value)

wref_mantel_table <- wref_eb_summary %>%
  left_join(wref_rad_select, by=c("TIMESTAMP", "group")) %>%
  mutate(RAD_dT = RAD_Tl - LAYER_TA) %>%
  select(-LAYER_TA, -RAD_Tl) %>%
  pivot_longer(contains("dT")) %>%
  pivot_wider(names_from="group", values_from="value")

## Mantel tests ----
old_wref_dT_dist <- old_wref_mantel_table %>%
  ungroup() %>%
  select(contains("canopy")) %>%
  dist()

old_wref_group_dist <- old_wref_mantel_table %>%
  ungroup() %>%
  mutate(is_TIR = temp_src == "TIR_dT") %>%
  pull(is_TIR) %>%
  dist()

# old_wref_mantel <- mantel(old_wref_dT_dist, old_wref_group_dist)
old_wref_mantel <- readRDS("data_out/502_analysis/old_wref_mantel")

wref_dT_dist <- wref_mantel_table %>%
  ungroup() %>%
  select(`0`:`3`) %>%
  dist()

wref_group_dist <- wref_mantel_table %>%
  ungroup() %>%
  mutate(is_rad = name == "RAD_dT") %>%
  pull(is_rad) %>%
  dist()

# wref_mantel <- mantel(wref_dT_dist, wref_group_dist)
wref_mantel <- readRDS("data_out/502_analysis/wref_mantel")

## Crossmatch tests ----

# old_wref_xmatch <- crossmatchtest(
#   old_wref_mantel_table$temp_src == "TIR_dT",
#   as.matrix(old_wref_dT_dist)
# )
# 
# wref_xmatch <- crossmatchtest(
#   wref_mantel_table$name == "RAD_dT",
#   as.matrix(wref_dT_dist)
# )
old_wref_xmatch <- readRDS("data_out/502_analysis/old_wref_xmatch")
wref_xmatch <- readRDS("data_out/502_analysis/wref_xmatch")

# Question 1: difference in dT as a function of envt + layer position ----
# Make response matrix
wref_eb_scale <- wref_eb %>%
  select(SITE_NEON, EB_MODEL_dT_LE, EB_MODEL_dT_Rn, EB_MODEL_omega,
         RAD_SW_ABS, RAD_LW_ABS, LAYER_TA, LAYER_RH, LAYER_WS, LAYER_L) %>%
  mutate(across(
    c(EB_MODEL_dT_LE, EB_MODEL_dT_Rn, EB_MODEL_omega, LAYER_L),
    function(x) (x - min(x)) / (max(x) - min(x))
  ))

resp_vars <- c("EB_MODEL_dT_LE", "EB_MODEL_dT_Rn", "EB_MODEL_omega")

wref_eb_dist <- wref_eb_scale %>%
  select(all_of(resp_vars)) %>%
  dist()

# Fit both models
# wref_perm_no_L <- adonis2(
#   wref_eb_dist ~ RAD_SW_ABS + RAD_LW_ABS + LAYER_TA + 
#     LAYER_RH + LAYER_WS,
#   data=wref_eb_scale
# )
# 
# wref_perm_L <- adonis2(
#   wref_eb_dist ~ RAD_SW_ABS + RAD_LW_ABS + LAYER_TA + 
#     LAYER_RH + LAYER_WS + LAYER_L,
#   data=wref_eb_scale
# )
wref_perm_no_L <- readRDS("data_out/502_analysis/wref_perm_no_L")
wref_perm_L <- readRDS("data_out/502_analysis/wref_perm_L")

# Calculate AICc

aicc_no_L <- AICc_permanova2(wref_perm_no_L)
aicc_L <- AICc_permanova2(wref_perm_L)
delta_aicc <- aicc_no_L[1, 1] - aicc_L[1, 1]

# Question 2: What drives more variance in T profiles, LE or Rg? ----
wref_eb_profiles <- wref_eb %>%
  select(TIMESTAMP, RAD_TOC_LW_IN, RAD_TOC_SW_IN, EB_MODEL_LE, 
         LAYER_L, LAYER_TA, EB_MODEL_Tl) %>%
  group_by(TIMESTAMP) %>%
  mutate(Rg = RAD_TOC_LW_IN + RAD_TOC_SW_IN,
         LE_total = sum(EB_MODEL_LE),
         dT = EB_MODEL_Tl - 273.15 - LAYER_TA) %>%
  ungroup() %>%
  select(TIMESTAMP, Rg, LE_total, dT, LAYER_L) %>%
  pivot_wider(names_from=LAYER_L, values_from=dT)

wref_eb_profile_dist <- wref_eb_profiles %>%
  select(`0`:`8`) %>%
  dist()

# Full data
wref_profile_perm <- adonis2(wref_eb_profile_dist ~ Rg + LE_total, 
                             data=wref_eb_profiles,
                             by="margin")

# Bootstrap
# set.seed(0221)
# n_iter <- 100
# wref_profile_perm_bootstrap <- foreach(i=1:n_iter, .combine=rbind) %do% {
#   wref_eb_profile_sample <- slice_sample(wref_eb_profiles, prop=0.9)
#   
#   wref_eb_profile_sample_dist <- wref_eb_profile_sample %>%
#     select(`0`:`8`) %>%
#     dist()
#   
#   wref_profile_sample_perm <- adonis2(wref_eb_profile_sample_dist ~ Rg + LE_total, 
#                                       data=wref_eb_profile_sample,
#                                       by="margin")
#   
#   c(
#     wref_profile_sample_perm[1, 3], wref_profile_sample_perm[2, 3]
#   )
# }
# 
# wref_profile_perm_boot_df <- wref_profile_perm_bootstrap %>%
#   as.data.frame() %>%
#   setNames(c("Rg", "LE"))

wref_profile_perm_boot_df <- read_csv("data_out/502_analysis/wref_profile_bootstrap.csv")

# Save results ----
# saveRDS(old_wref_mantel, "data_out/502_analysis/old_wref_mantel")
# saveRDS(wref_mantel, "data_out/502_analysis/wref_mantel")
# saveRDS(old_wref_xmatch, "data_out/502_analysis/old_wref_xmatch")
# saveRDS(wref_xmatch, "data_out/502_analysis/wref_xmatch")
# saveRDS(wref_perm_L, "data_out/502_analysis/wref_perm_L")
# saveRDS(wref_perm_no_L, "data_out/502_analysis/wref_perm_no_L")
# saveRDS(wref_profile_perm, "data_out/502_analysis/wref_profile_perm")
# write_csv(wref_profile_perm_boot_df, "data_out/502_analysis/wref_profile_bootstrap.csv")

# Figures ----
pres_theme <- theme(
  axis.text = element_text(size=20),
  strip.text = element_text(size=20),
  axis.title = element_text(size=20)
)

## 1:1 plots ----
old_wref_1to1 <- old_wref_model_tir %>%
  ungroup() %>%
  mutate(REF_DT = TIR_Tl - LAYER_TA,
         REF_TL = TIR_Tl,
         REF="Thermal camera",
         EB_TL = EB_MODEL_Tl,
         EB_DT = EB_MODEL_Tl - LAYER_TA) %>%
  select(group, LAYER_TA, EB_TL, EB_DT, REF, REF_TL, REF_DT)

wref_1to1 <- wref_eb_summary %>%
  left_join(wref_rad_select, by=c("TIMESTAMP", "group")) %>%
  ungroup() %>%
  mutate(RAD_dT = RAD_Tl - LAYER_TA,
         group = case_match(
           group,
           0 ~ "Upper canopy",
           1 ~ "Middle canopy",
           2 ~ "Middle canopy",
           3 ~ "Lower canopy"
         ),
         REF_DT = RAD_dT,
         REF_TL = RAD_Tl,
         REF="Radiometer",
         EB_TL=EB_MODEL_dT + LAYER_TA,
         EB_DT=EB_MODEL_dT
  ) %>%
  select(group, LAYER_TA, EB_TL, EB_DT, REF, REF_TL, REF_DT)

old_wref_1to1 %>%
  bind_rows(wref_1to1) %>%
  mutate(group = factor(group, levels=c("Upper canopy", "Middle canopy", "Lower canopy"))) %>%
  ggplot(aes(x=REF_DT, y=EB_DT)) +
  geom_point(alpha=0.5, stroke=NA) +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  facet_grid(group ~ REF) +
  theme_bw() +
  labs(x=expression("Reference"~Delta*"T (K)"),
       y=expression("Model"~Delta*"T (K)"))

# Same but leaf temperature instead of delta T
old_wref_1to1 %>%
  bind_rows(wref_1to1) %>%
  mutate(group = factor(group, levels=c("Upper canopy", "Middle canopy", "Lower canopy"))) %>%
  ggplot(aes(x=REF_TL, y=EB_TL)) +
  geom_point(alpha=0.5, stroke=NA) +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  facet_grid(group ~ REF) +
  theme_bw() +
  labs(x=expression("Reference T"[L]~"("*degree*"C)"),
       y=expression("Model T"[L]~"("*degree*"C)")) +
  pres_theme

## Vertical gradients ----
wref_eb %>%
  select(LAYER_L, PS_LAYER_GPP, EB_MODEL_LE, EB_MODEL_omega) %>%
  mutate(LAYER_L = factor(LAYER_L)) %>%
  pivot_longer(-LAYER_L) %>%
  mutate(name = case_match(name,
      "EB_MODEL_LE" ~ "Latent heat flux (W m-2)",
      "EB_MODEL_omega" ~ "Decoupling coefficient",
      "PS_LAYER_GPP" ~ "Photosynthesis (umol m-2 s-1)"
  )) %>%
  ggplot(aes(x=factor(LAYER_L), y=value)) +
  geom_boxplot(aes(fill=name)) +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  facet_wrap(~ name, strip.position="bottom", scales="free_x") +
  theme_bw() +
  labs(y="", x="") +
       #x=expression("Cumulative LAI"~"(m"^2~"m"^-2*")")) +
  theme(legend.position="none",
        strip.placement="outside",
        strip.background = element_rect(fill=NA, linewidth=NA),
        strip.text = element_text(size=10),
        panel.spacing = unit(2, "lines")) +
  pres_theme +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_blank())

## Variation explained histograms ----
wref_profile_perm_boot_df %>%
  pivot_longer(everything()) %>%
  mutate(name = case_match(name,
                          "LE" ~ "Transpiration",
                          "Rg" ~ "Downwelling radiation"
                          )) %>%
  ggplot(aes(x=value)) +
  geom_histogram(aes(fill=name)) +
  labs(y="Count",
       x="Type III proportion of variance explained",
       fill="") +
  theme_bw() +
  theme(legend.position=c(1, 0.9),
        legend.justification = "right",
        legend.background = element_rect(fill=NA)) +
  pres_theme +
  theme(legend.text = element_text(size=14))
  
## Plot of model II slopes by layer ----
wref_eb_linmods <- wref_eb %>%
  select(LAYER_L, LAYER_TA, EB_MODEL_Tl) %>%
  mutate(EB_MODEL_Tl = EB_MODEL_Tl - 273.15) %>%
  group_by(LAYER_L) %>%
  nest() %>%
  mutate(
    mod2 = map(data, function(d) lmodel2(EB_MODEL_Tl ~ LAYER_TA, data=d)),
    slope_p025 = map_dbl(mod2, function(m) m$confidence.intervals[3, 4]),
    slope_p50  = map_dbl(mod2, function(m) m$regression.results[3, 3]),
    slope_p975 = map_dbl(mod2, function(m) m$confidence.intervals[3, 5]),
    int_p50    = map_dbl(mod2, function(m) m$regression.results[3, 2])
  ) 

wref_eb_linmods %>%
  ggplot(aes(x=slope_p50, y=LAYER_L)) +
  geom_vline(xintercept=1, color="grey50", linetype="dashed") +
  geom_errorbar(aes(y=LAYER_L, xmin=slope_p025, xmax=slope_p975)) +
  geom_point(aes(fill=int_p50), color="black", size=3, pch=21) +
  scale_y_reverse(breaks=0:9) +
  scale_fill_viridis_c() +
  labs(x=expression("T"[leaf]/"T"[air]~"model II slope"~"("*degree*"C/"*degree*"C)"),
       y=expression("Cumulative leaf area index (m"^2~"m"^-2*")"),
       fill=expression("Model II intercept ("*degree*"C)")) +
  theme_bw() +
  guides(fill = guide_colorbar(
    ticks.colour="black",
    frame.colour="black"
  )) +
  pres_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))

## Proportion shaded GPP calculation ----
wref_eb %>%
  mutate(
    LIDAR_prop_sunlit = exp(LIDAR_kb_black * LAYER_L),
    PS_SUNLIT_GPP = LIDAR_prop_sunlit * PS_LAYER_GPP,
    PS_SHADE_GPP = (1-LIDAR_prop_sunlit) * PS_LAYER_GPP
  ) %>%
  summarize(SUNLIT_GPP = sum(PS_SUNLIT_GPP),
            SHADE_GPP = sum(PS_SHADE_GPP))

wref_eb_cool <- wref_eb %>%
  filter(EB_MODEL_Tl - 273.15 - LAYER_TA < 0)

## Density plots of mantel distances
distance_df <- data.frame(
  temp_dist=as.numeric(wref_dT_dist),
  group_dist=as.numeric(wref_group_dist),
  period="Radiometer"
) %>% bind_rows(data.frame(
  temp_dist=as.numeric(old_wref_dT_dist),
  group_dist=as.numeric(old_wref_group_dist),
  period="Thermal camera"
))

distance_df %>%
  mutate(
    group_dist = case_match(
      group_dist,
      0 ~ "Model vs. Model &\nSensor vs. Sensor",
      1 ~ "Model vs. Sensor"
    )
  ) %>%
  ggplot(aes(x=temp_dist, y=group_dist)) +
  geom_density_ridges(aes(fill=group_dist)) +
  facet_wrap(~ period, scales="free_x") +
  labs(x="Temperature profile Euclidean distance (K)",
       y="") +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_text(size=20),
        strip.text=element_text(size=20),
        axis.text=element_text(size=20))
