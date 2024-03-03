library(ggplot2)
library(ggrepel)
library(latex2exp)
library(tidyverse)

theme_set(ggdark::dark_theme_bw())

# Resolution tradeoff figure ----
tsats <- data.frame(
  name=c("MODIS", "Sentinel-3", "ECOSTRESS", "Landsat", "NASA G-LiHT", "NEON AOP"),
  spatial=c(1000, 1000, 70, 30, 0.1, 0.1),
  temporal=c(24, 24, 24, 192, 8760, 8760)
)

anno_color <- "lightblue"
anno_size <- 12

ggplot(tsats, aes(x=temporal, y=spatial)) +
  geom_point() +

  geom_smooth(method="lm", se=FALSE, linetype="dashed", color="gray") +
  annotate("text", x=1000, y=4, label="← Resolution Tradeoff →", angle=-40,
           color="gray", size=anno_size) +
  # Annotations for LiDAR and flux tower resolution
  # geom_vline(xintercept=0.5, color=anno_color) +
  # geom_hline(yintercept=0.1, color=anno_color) +
  # annotate("text", y=0.15, x=10, label="LiDAR",
  #          color=anno_color, size=anno_size, fontface="italic", alpha=0.5) +
  # annotate("text", x=0.55, y=10, label="Flux tower",
  #          color=anno_color, hjust="left", size=anno_size, fontface="italic",
  #          alpha=0.5) +
  # annotate("point", x=0.5, y=0.1, pch=10, color=anno_color, size=4) +
  geom_label_repel(aes(label=name), size=anno_size - 2) +
  labs(x="Temporal Resolution",
       y="Spatial Resolution") +
  scale_x_log10(breaks=c(2.4, 24, 240, 2400),
                labels=c("Hourly", "Daily", "Weekly", "Yearly"),
                limits=c(10^-0.31, 10^4)) +
  scale_y_log10(breaks=10^(-1:3),
                labels=sapply(-1:3, function(exp) TeX(paste0("10$^{", exp, "}$ m")))) +
  #ggtitle("Our model overcomes the resolution tradeoff in traditional thermal data") +
  theme(
    axis.text = element_text(size=36),
    axis.title = element_text(size=40)
  )

#ggsave("graphics/resolution_figure2.png", width=19.2, height=10.8, units="in",
#       dpi=600)

ggplot(tsats, aes(x=temporal, y=spatial)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", color="gray") +
  geom_label_repel(aes(label=name)) +
  annotate("text", x=1000, y=4, label="← Resolution Tradeoff →", angle=-37,
           color="gray", size=anno_size) +
  labs(x="Temporal Resolution",
       y="Spatial Resolution") +
  scale_x_log10(breaks=c(2.4, 24, 240, 2400),
                labels=c("Hourly", "Daily", "Weekly", "Yearly"),
                limits=10^c(-0.31, 4)) +
  scale_y_log10(breaks=10^(-1:3),
                labels=sapply(-1:3, function(exp) TeX(paste0("10$^{", exp, "}$ m"))),
                limits=10^c(-1, 3)) +
  #ggtitle("Our model overcomes the resolution tradeoff in traditional thermal data") +
  theme(
    axis.text = element_text(size=18),
    axis.title = element_text(size=20)
  )

# Freilich figures ----

library(tidyverse)
library(ggdark)

theme_set(dark_theme_bw())

big_theme <- theme(
  axis.text = element_text(size=28),
  axis.title = element_text(size=32),
  legend.title = element_text(size=28),
  legend.text = element_text(size=28),
  strip.text = element_text(size=24),
  panel.grid = element_line(linewidth=1, color=rgb(.2, .2, .2))
)

all_tc_result <- read_csv("data_out/cross_site_gs_tc_result.csv") %>%
  filter(z > 1) %>%
  group_by(SITE) %>%
  mutate(
    position = case_when(
      abs(z - max(z)) < 1 ~ "Crown",
      abs(z - min(z)) < 1 ~ "Lower",
      .default="Middle"
    ),
    position = factor(position, c("Crown", "Middle", "Lower"))
  ) %>% ungroup()


## 1:1 plot with bias annotation ----
bias_df <- all_tc_result %>%
  group_by(SITE, position) %>%
  summarize(bias = mean(t_canopy_model_still_Tl - t_canopy, na.rm=TRUE)) %>%
  mutate(bias_pretty = str_c("Bias = ", format(round(bias, digits=2)), " °C"))

all_tc_result %>%
  #filter(SITE %in% c("US-xWR", "US-xTE")) %>%
  #filter(t_canopy < (T_lw - 273.15) + 10) %>%
  group_by(SITE) %>%
  #sample_frac(0.5) %>%
  ggplot(aes(x=t_canopy, y=t_canopy_model_still_Tl)) +
  geom_bin2d(bins=50) +
  geom_abline(slope=1, intercept=0) +
  # geom_text(
  #   data=bias_df %>% filter(SITE %in% c("US-xTE", "US-xWR")), 
  #   mapping=aes(label=bias_pretty),
  #   x=15, y=40, size=6
  # ) +
  coord_equal() +
  scale_fill_distiller(palette="RdYlBu", limits=c(0, 250), oob=scales::squish) +
  facet_grid(position ~ SITE) +
  labs(x="Tower radiometer temperature (°C)",
       y="Modeled leaf temperature (°C)",
       fill="Count") +
  big_theme +
  theme(strip.text = element_text(size=18),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))

#ggsave("graphics/poster_one_to_one.png", dpi=300, height=10.8)

## Mean LE/Rn forcing ----
all_tc_forcing <- all_tc_result %>%
  mutate(
    Rn_forcing = t_notransp_model_still_Tl - ta_interp,
    LE_forcing = t_canopy_model_still_Tl - t_notransp_model_still_Tl,
    dT = Rn_forcing + LE_forcing
  )

forcing_means <- all_tc_forcing %>%
  select(SITE, position, Rn_forcing, LE_forcing, dT) %>%
  pivot_longer(-c(SITE, position)) %>%
  group_by(SITE, position, name) %>%
  summarize(median = median(value, na.rm=TRUE),
            sdev = sd(value, na.rm=TRUE))

forcing_means %>%
  filter(name != "dT") %>%
  ggplot(aes(x=position, y=median)) +
  geom_bar(aes(fill=name), stat="identity") +
  #geom_errorbar(aes(color=name, ymin=median - sdev, ymax=median + sdev)) +
  geom_point(
    data=forcing_means %>% filter(name == "dT"),
    mapping=aes(y=median, x=position, color="Net"),
    size=4
  ) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  scale_fill_manual(
    values=c(
    "LE_forcing"="#158FFF",
    "Rn_forcing"="#FF1919"
    ),
    labels=c(expression(paste(lambda, "E")), "Rn")
  ) +
  scale_color_manual(
    values=c("Net"="white")
  ) +
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 2, title=NULL)) +
  facet_wrap(~ SITE) +
  labs(
    x=expression(paste(T[leaf], " - ", T[air], " (K)")), 
    y="", fill="", color=""
  )# +
  big_theme

#ggsave("graphics/freilich_mean_forcing.png", dpi=300, height=10.8, width=16)

## Second version to explain the plot ----
fake_forcing <- data.frame(
  name=c("Rn_forcing", "LE_forcing"),
  value=c(1, -2)
)

ggplot(fake_forcing, aes(y=value, fill=name)) + 
  geom_bar(aes(x=1), stat="identity") +
  coord_flip()

fake_forcing %>%
  ggplot(aes(x=value, y=1)) +
  geom_bar(aes(fill=name), stat="identity") + 
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  scale_fill_manual(
    values=c(
      "LE_forcing"="#158FFF",
      "Rn_forcing"="#FF1919"
    ),
    labels=c(expression(paste(lambda, "E")), "Rn")
  ) +
  scale_color_manual(
    values=c("Net"="white")
  ) +
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 2, title=NULL)) +
  labs(
    x=expression(paste(T[leaf], " - ", T[air], " (K)")), 
    y="", fill="", color=""
  )
  big_theme

## Proportion of variance explained ----
all_tc_linmods <- all_tc_result %>%
  mutate(
    Rn_forcing = t_notransp_model_still_Tl - ta_interp,
    LE_forcing = t_canopy_model_still_Tl - t_notransp_model_still_Tl,
    dT = Rn_forcing + LE_forcing
  ) %>%
  group_by(SITE) %>%
  select(SITE, t_canopy_model_still_Tl, dT, Rn_forcing, LE_forcing, ta_interp, 
         ws_interp, vpd, sw_in, gs_interp) %>%
  nest(data=-c(SITE)) %>%
  mutate(
    anova_forcing = map(data, function (d) {
      anova(lm(
        t_canopy_model_still_Tl ~ ta_interp + Rn_forcing + LE_forcing, data=d
      ))
    }),
    anova_dt = map(data, function(d) {
      anova(lm(dT ~ ta_interp + ws_interp + vpd + gs_interp + sw_in, data=d))
    }),
    varexp_forcing = map(anova_forcing, function(av) {
      avss <- av$"Sum Sq"
      prop_exp <- avss / sum(avss)
      
      data.frame(var = rownames(av), prop_exp=prop_exp)
    }),
    varexp_dt = map(anova_dt, function(av) {
      avss <- av$"Sum Sq"
      prop_exp <- avss / sum(avss)
      
      data.frame(var = rownames(av), prop_exp=prop_exp)
    })
  )

# ... in T_leaf
all_tc_linmods %>%
  select(SITE, varexp_forcing) %>%
  unnest(varexp_forcing, names_sep="_") %>%
  filter(varexp_forcing_var != "Residuals") %>%
  ggplot(aes(x=SITE, y=varexp_forcing_prop_exp)) +
  geom_bar(aes(fill=varexp_forcing_var), stat="identity", position="dodge") +
  facet_wrap(~ SITE, nrow=3, scales="free_x") +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="Proportion of variance", fill="") +
  scale_fill_manual(
    values=c(
      "LE_forcing"="#158FFF",
      "Rn_forcing"="#FF1919",
      "ta_interp"="#DDDDDD"
    ),
    labels=c(expression(paste(lambda, "E")), "Rn", "Air temperature")
  ) +
  big_theme +
  theme(legend.text.align = 0)

#ggsave("graphics/freilich_tleaf_var.png", dpi=300, height=10.8, width=16)

# ... in dT
all_tc_linmods %>%
  select(SITE, varexp_dt) %>%
  unnest(varexp_dt, names_sep="_") %>%
  filter(varexp_dt_var != "Residuals") %>%
  mutate(
    varexp_dt_var = case_match(
      varexp_dt_var,
      "gs_interp" ~ "Stomatal conductance",
      "sw_in" ~ "Solar radiation",
      "ta_interp" ~ "Air temperature",
      "vpd" ~ "Vapor pressure deficit",
      "ws_interp" ~ "Wind speed"
    )
  ) %>%
  ggplot(aes(x=SITE, y=varexp_dt_prop_exp)) +
  geom_bar(aes(fill=varexp_dt_var), stat="identity", position="dodge") +
  facet_wrap(~ SITE, nrow=3, scales="free_x") +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="Proportion of variance", fill="") +
  ylim(0, 1) +
  big_theme +
  theme(legend.text.align = 0)

#ggsave("graphics/freilich_dt_var.png", dpi=300, height=10.8, width=16)


# Site map ----
library(maps)
library(mapdata)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(sf)
library(sp)

theme_set(theme_bw())

site_info <- read_csv("data_working/neon_site_metadata.csv") %>%
  filter(!(site_neon %in% c("YELL", "SOAP"))) %>%
  st_as_sf(coords=c("tower_lon", "tower_lat"), crs=st_crs(4326)) %>%
  cbind(., st_coordinates(.))

sites_l48 <- site_info %>% filter(Y < 49)
sites_ak  <- site_info %>% filter(Y > 49)
stopifnot(nrow(sites_l48) + nrow(sites_ak) == nrow(site_info))

# N.b. it is challenging to plot lower 48 and Alaska in their true position,
# so we use insets.
states <- read_sf("data_in/cb_2022_us_state_20m/cb_2022_us_state_20m.shp")

ak  <- states %>% filter(STUSPS == "AK")
hi  <- states %>% filter(STUSPS == "HI")
# sorry, Puerto Rico :(
l48 <- states %>% filter(!(STUSPS %in% c("AK", "HI", "PR")))

label_size <- 10

# Inset maps
(
  ak_map <- ggplot(data=ak) +
    geom_sf() +
    geom_label_repel(
      data=sites_ak, 
      mapping=aes(label=site_ameriflux, geometry=geometry),
      stat="sf_coordinates",
      size=label_size
    ) +
    coord_sf(crs=st_crs(3467)) +
    theme_void()
)

# Now the mainland map
(
  mainland <- ggplot(data=l48) +
    geom_sf() +
    geom_label_repel(
      data=sites_l48, 
      mapping=aes(label=site_ameriflux, geometry=geometry),
      stat="sf_coordinates",
      size=label_size
    ) +
    coord_sf(crs=st_crs(9311), 
             xlim = c(-2300000, 2600000), ylim = c(-2500000, 740000),
             expand=FALSE, datum=NA) +
    labs(x="", y="")
)

# Add them all together
ak_xmin <- -2200000
ak_ymin <- -2400000
scale <- 2.5

ak_xmax <- ak_xmin + (1600000 - (-2400000))/scale
ak_ymax <- ak_ymin + (2500000 - 200000)/scale

mainland +
  annotation_custom(
    grob = ggplotGrob(ak_map),
    xmin = ak_xmin,
    xmax = ak_xmin + (1600000 - (-2400000))/scale,
    ymin = ak_ymin,
    ymax = ak_ymin + (2500000 - 200000)/scale
  ) +
  annotate(
    "line",
    x=c(ak_xmin, ak_xmax-150000, ak_xmax+200000),
    y=c(ak_ymax-50000, ak_ymax-50000, ak_ymin),
    color="white"
  ) +
  labs(x="", y="") +
  theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_)
  )
  
# ggsave("graphics/site_map.png", dpi=300, width=19.20, height=10.80,
#        bg="transparent")


# Tower/canopy heights ----
library(ggpattern)
library(ggdark)

tower_color="grey50"

site_info <- site_info %>%
  mutate(filepath = file.path(
    "graphics", "tree_svgs", str_c(str_to_lower(site_ameriflux), ".png")
  ))

radiometer_positions <- all_tc_result %>%
  select(SITE, z) %>% distinct()

file_match_vector <- site_info$filepath
names(file_match_vector) <- site_info$site_ameriflux

ggplot(site_info, aes(x=site_ameriflux, y=tower_height)) +
  # Add the tower top first, this can make it easier to not clip
  # the tower top if the image is resized.
  geom_bar_pattern(
    aes(y=tower_height+2),
    fill="black",
    width=0.3,
    stat="identity",
    pattern="image",
    pattern_type = "none",
    pattern_filename="graphics/tower_top_light.png",
    pattern_scale=-1,
    pattern_gravity="north",
    pattern_yoffset=0.5,
    alpha=0
  ) +
  # Tower tile
  geom_bar_pattern(
    aes(y=tower_height-5),
    stat="identity",
    fill="black",
    pattern="image",
    pattern_type="tile",
    pattern_filename="graphics/tower_tile_light.png",
    pattern_filter="box",
    pattern_scale=-1,
    color=tower_color,
    linewidth=1,
    width=0.2,
    alpha=0
  ) +
  # Tree outlines
  # Now add the trees
  geom_bar_pattern(
    aes(pattern_filename=site_ameriflux, y=canopy_height),
    stat="identity",
    pattern="image",
    pattern_gravity="south",
    alpha=0
  ) +
  scale_pattern_filename_manual(values=file_match_vector) +
  # Debug geom to make sure trees match the actual height
  # geom_point(aes(y=canopy_height), color="green") +
  # Radiometer heights
  # geom_point(
  #   mapping=aes(x=SITE, y=z), data=radiometer_positions,
  #   color="red"
  # ) +
  # Suppress image legend
  guides(pattern_filename="none") +
  labs(x="", y="Height above ground (m)") +
  dark_theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="white"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=16))

# ggsave("graphics/tower_canopy_height.png", width=9.5, height=4)

# WREF sensor positions ----
new_var_heights <- amerifluxr::amf_var_info() %>%
  filter(Site_ID == "US-xWR",
         str_detect(Variable, 
                    "H2O_\\d_\\d_\\d|PPFD_IN_\\d_\\d_\\d|TA_\\d_\\d_\\d|WS_\\d_\\d_\\d")) %>%
  mutate(var_type = str_split_i(Variable, "_", 1),
         var_type = case_match(
           var_type,
           "H2O" ~ "Humidity",
           "PPFD" ~ "Light",
           "TA" ~ "Air temp.",
           "WS" ~ "Wind speed"
         ),
         Height_floor = floor(Height/2)*2)

old_var_heights <- data.frame(
  var_type=rep(c("Air temp.", "Humidity", "Wind speed", "Light"), 2),
  Height=c(rep(2, 4), rep(70, 4))
)
  
wref_zh <- 58
wref_zr <- 74

ggplot() +
  # Tree outline
  # geom_bar_pattern(
  #   aes(x=1, y=wref_zh),
  #   pattern_filename="graphics/tree_svgs/us-xwr.png",
  #   stat="identity",
  #   pattern="image",
  #   alpha=0
  # ) +
  # # Tower top
  # geom_bar_pattern(
  #   aes(x=2, y=wref_zr),
  #   fill="black",
  #   width=0.2,
  #   stat="identity",
  #   pattern="image",
  #   pattern_type = "none",
  #   pattern_filename="graphics/tower_top.png",
  #   pattern_scale=-1,
  #   pattern_gravity="north",
  #   pattern_yoffset=0.5,
  #   alpha=0
  # ) +
  # # Tower tile
  # geom_bar_pattern(
  #   aes(x=2, y=wref_zr-8),
  #   stat="identity",
  #   fill="black",
  #   color="black",
  #   pattern="image",
  #   pattern_type="tile",
  #   pattern_filename="graphics/tower_tile.png",
  #   pattern_filter="box",
  #   pattern_scale=-1,
  #   linewidth=1,
  #   width=0.1,
  #   alpha=0
  # ) +
  # Horizontal dashed line
  geom_hline(yintercept=wref_zh, linetype="dashed", color="grey50") +
  annotate("text", x=-Inf, y=wref_zh, vjust=-1, hjust=0, label=" Canopy height",
           fontface="italic", color="grey50") +
  # Modern sensor positions
  geom_point(
    mapping=aes(x="2018-Present", y=Height, color=var_type),
    data=new_var_heights,
    position=position_dodge(width=0.1)
  ) +
  # Old sensor positions
  geom_point(
    mapping=aes(x="2014-2015", y=Height, color=var_type),
    data=old_var_heights,
    position=position_dodge(width = 0.1)
  ) +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank()) +
  labs(x="", y="Height above ground (m)",
       color="Sensor type")
