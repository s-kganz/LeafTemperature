library(lidR)
library(tidyverse)

site_lai <- read_csv("data_out/neon_lai.csv") %>%
  filter(!(site %in% c("YELL", "SOAP"))) %>%
  group_by(site) %>%
  summarize(LAI = median(L, na.rm=TRUE))

site_info <- read_csv("data_working/neon_site_metadata.csv") %>%
  left_join(site_lai, by=c("site_neon"="site")) %>%
  drop_na()

tc_data <- read_csv("data_out/cross_site_tc_data.csv")

fit_lidar_constants <- function(site, site_lai, zmax=NA, zmin=1) {
  # Z where cumulative LAI == total LAI
  # Measurement height of Io
  if (is.na(zmax)) {
    zmax <- tc_data %>% filter(SITE_NEON == site) %>% 
      pull(z) %>% max() %>% ceiling()
  }
  
  las <- readLAS(file.path(
    "data_working/neon_neartower_lidar", str_c(site, ".laz")
  )) %>%
    filter_poi(Z > zmin, Z < zmax)
  
  # Calibrate k so that the integral of LAD == LAI. Reverse the order of rows
  # so we can compute the cumulative LAI looking *down* the canopy.
  lad_uncal_sum <- LAD(las$Z, k=1, z0 = zmin) %>% pull(lad) %>% sum(na.rm=TRUE)
  k_lai <- lad_uncal_sum / site_lai
  
  lad_profile <- LAD(las$Z, k=k_lai, z0 = zmin)
  lad_profile <- lad_profile[nrow(lad_profile):1, ]
  lad_profile$cum_lai <- cumsum(lad_profile$lad)
  
  # Now calculate transmission profile through the canopy. For each z, calc geometric
  # mean transmittance to that point in the canopy.
  tc_site <- tc_data %>% filter(SITE_NEON == site)
  
  iv_io_df <- tc_site %>%
    # Drop night time data
     filter(SW_IN_TOC > 100) %>%
    select(TIMESTAMP, PPFD_IN_TOC, z, zmax, par_interp) %>%
    group_by(TIMESTAMP) %>%
    mutate(iv = par_interp,
           io = PPFD_IN_TOC) %>%
    ungroup() %>%
    mutate(iv_io = iv / io) %>%
    group_by(z) %>%
    summarize(
      iv_io_mean = suppressWarnings(log(iv_io)) %>% mean(na.rm=TRUE) %>% exp(),
      # https://en.wikipedia.org/wiki/Geometric_standard_deviation
      iv_io_sd = suppressWarnings(log(iv_io / iv_io_mean)^2) %>% 
        mean(na.rm=TRUE) %>% sqrt() %>% exp()
    ) %>%
    mutate(cum_lai = approx(lad_profile$z, lad_profile$cum_lai, z, rule=2)$y,
           site = site)
  
  # Force the intercept at 0
  iv_io_lm <- lm(I(log(iv_io_mean) - 1) ~ cum_lai + 0, data=iv_io_df)
  
  return(list(
    data=iv_io_df,
    constants=list(
      site = site,
      k_lai = k_lai,
      k_rad = -iv_io_lm$coefficients[1]
    )
  ))
}

lidar_constants <- lapply(
  1:nrow(site_info),
  function(i) {
    fit_lidar_constants(
      site=site_info$site_neon[i],
      site_lai=site_info$LAI[i],
      zmax=site_info$canopy_height[i],
      zmin=1
    )
  }
)

iv_io_data <- lapply(lidar_constants, function(l) l$data) %>% bind_rows()
lidar_constants_df <- lapply(lidar_constants, function(l) l$constants) %>% bind_rows()

# Plot mean transmission points over the fitted line
iv_io_lines <- expand.grid(
    cum_lai=seq(0, 5, length.out=20), 
    site=lidar_constants_df$site
  ) %>%
  mutate(
    k_rad = lidar_constants_df$k_rad[match(site, lidar_constants_df$site)],
    iv_io = exp(-k_rad * cum_lai)
  )

iv_io_data %>%
  ggplot(aes(x=cum_lai, y=iv_io_mean)) +
  geom_point() +
  geom_line(data=iv_io_lines, aes(x=cum_lai, y=iv_io), color="red") +
  facet_wrap(~ site) +
  scale_y_log10() +
  theme_bw()

#write_csv(lidar_constants_df, "data_out/cross_site_lidar_constants.csv")

osbs_ivio <- tc_data %>%
  filter(SITE_NEON == "OSBS") %>%
  select(TIMESTAMP, SW_IN_TOC, z, PPFD_IN_TOC, par_interp) %>% 
  group_by(TIMESTAMP) %>% 
  mutate(iv = par_interp, io = PPFD_IN_TOC) %>% 
  ungroup() %>% mutate(iv_io = iv / io) %>%
  mutate(TIMESTAMP_LOCAL = TIMESTAMP + hours(-5))

osbs_ivio %>%
  mutate(hour = hour(TIMESTAMP_LOCAL)) %>%
  ggplot(aes(x=hour, y=SW_IN_TOC)) +
  geom_boxplot(aes(group=hour))

osbs_ivio %>%
  ggplot(aes(x=iv_io, y=z)) +
  geom_boxplot(aes(group=z)) +
  geom_vline(xintercept=1)

  