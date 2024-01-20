# To run the radiation submodel, we need empirical light attenutation
# coefficients for direct and diffuse radiation. These are calculated empirically
# based on observations where diffuse and direct light dominate.

library(lidR)
library(tidyverse)

site_lai <- read_csv("data_out/neon_sampled_dhp_lai.csv")

# If there is a prior estimate of LAI in the literature, use it. Otherwise,
# use the DHP-derived estimate.
site_info <- read_csv("data_working/neon_site_metadata.csv") %>%
  left_join(site_lai, by=c("site_neon"="site")) %>%
  mutate(
    lai_best = ifelse(is.na(lai), L_dhp_corrected, lai)
  )

tc_data <- read_csv("data_out/cross_site_tc_data.csv") %>%
  mutate(proportion_diffuse = SW_DIF / (SW_IN_TOC + SW_DIF),
         proportion_diffuse = pmin(pmax(proportion_diffuse, 0), 1))

fit_lidar_constants <- function(site, site_lai, zmax=NA, zmin=1) {
  # Z where cumulative LAI == total LAI and measurement height of Io
  if (is.na(zmax)) {
    zmax <- tc_data %>% filter(SITE_NEON == site) %>% 
      pull(z) %>% max() %>% ceiling()
  }
  
  las <- readLAS(file.path(
    "data_working/neon_neartower_lidar", str_c(site, ".laz")
  )) %>%
    filter_poi(Z > zmin, Z < zmax)
  
  # Calibrate k so that the integral of LAD == LAI. 
  lad_uncal_sum <- LAD(las$Z, k=1, z0 = zmin) %>% pull(lad) %>% sum(na.rm=TRUE)
  k_lai <- lad_uncal_sum / site_lai
  
  lad_profile <- LAD(las$Z, k=k_lai, z0 = zmin)
  # Reverse the order of rows so we can compute the cumulative LAI looking 
  # *down* the canopy.
  lad_profile <- lad_profile[nrow(lad_profile):1, ]
  lad_profile$cum_lai <- cumsum(lad_profile$lad)
  
  # Now calculate transmission profile through the canopy under direct and
  # diffuse conditions. For each z, calc geometric mean transmittance to that 
  # point in the canopy.
  tc_site <- tc_data %>% filter(SITE_NEON == site, SW_IN_TOC > 100) %>%
    select(TIMESTAMP, PPFD_IN_TOC, z, zmax, par_interp, proportion_diffuse)
  
  iv_io_df <- tc_site %>%
    mutate(iv = par_interp,
           io = PPFD_IN_TOC,
           proportion_diffuse_bin = round(proportion_diffuse, digits=1)) %>%
    mutate(iv_io = iv / io) %>%
    group_by(z, proportion_diffuse_bin) %>%
    summarize(
      iv_io_mean = suppressWarnings(log(iv_io)) %>% median(na.rm=TRUE) %>% exp(),
      # https://en.wikipedia.org/wiki/Geometric_standard_deviation
      iv_io_sd = suppressWarnings(log(iv_io / iv_io_mean)^2) %>% 
        mean(na.rm=TRUE) %>% sqrt() %>% exp()
    ) %>%
    drop_na() %>%
    mutate(cum_lai = approx(lad_profile$z, lad_profile$cum_lai, z, rule=2)$y,
           site = site)
  
  # Force the y-intercept at log(iv_io) = 0
  iv_io_lm <- lm(
    log(iv_io_mean) ~ cum_lai + cum_lai:proportion_diffuse_bin + 0, 
    data=iv_io_df
  )
  
  return(list(
    site=site,
    data=iv_io_df,
    rad_coef=coef(iv_io_lm),
    lai_coef=k_lai
  ))
}

lidar_constants <- lapply(
  1:nrow(site_info),
  function(i) {
    fit_lidar_constants(
      site=site_info$site_neon[i],
      site_lai=site_info$lai_best[i],
      zmax=site_info$canopy_height[i],
      zmin=1
    )
  }
)

iv_io_data <- lapply(lidar_constants, function(l) l$data) %>% bind_rows()
lidar_constants_df <- lapply(lidar_constants, function(l) {
  list(
    site=l$site,
    lai_coef=l$lai_coef,
    # Force attenuation coefficients less than 0
    kb_real=min(unname(l$rad_coef[1]), 0),
    kd_real=min(unname(l$rad_coef[1] + l$rad_coef[2]), 0)
  )
}) %>% bind_rows()

# Extrapolate lines for beam/diffuse coefficients
clai <- seq(0, 7, length.out=20)
iv_io_extrapolate <- expand.grid(clai, 1:nrow(lidar_constants_df)) %>%
  cbind(lidar_constants_df[.$Var2, ]) %>%
  select(-Var2) %>% rename(clai = Var1) %>%
  mutate(iv_io_beam = exp(kb_real * clai),
         iv_io_diff = exp(kd_real * clai)) %>%
  pivot_longer(contains("iv_io"))

# Plot mean transmission points over the fitted line
ggplot(iv_io_data, aes(x=cum_lai, y=iv_io_mean)) +
  geom_point(aes(color=proportion_diffuse_bin)) +
  geom_line(data=iv_io_extrapolate, mapping=aes(x=clai, y=value, group=name, linetype=name)) +
  facet_wrap(~ site, scales="free") +
  labs(x="Cumulative LAI", y="Proportion available light",
       color="Proportion diffuse light", name="") +
  scale_y_log10()

write_csv(lidar_constants_df, "data_out/neon_lidar_constants.csv")

  