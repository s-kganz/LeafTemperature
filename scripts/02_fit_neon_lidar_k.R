# To run the radiation submodel, we need empirical light attenutation
# coefficients for direct and diffuse radiation. These are calculated empirically
# based on observations where diffuse and direct light dominate. We also need
# to map vertical heights above ground to cumulative leaf area index. This
# script calculates the attenuation coefficients and generates LAI-Z crosswalk
# tables.

library(lidR)
library(tidyverse)
library(amerifluxr)

source("scripts/tower_util.R")

# Prepare site-level LAIs ----
site_lai <- read_csv("data_out/neon_sampled_dhp_lai.csv")

# If there is a prior estimate of LAI in the literature, use it. Otherwise,
# use the DHP-derived estimate.
site_info <- read_csv("data_working/neon_site_metadata.csv") %>%
  left_join(site_lai, by=c("site_neon"="site")) %>%
  mutate(
    lai_best = ifelse(is.na(lai), L_dhp_corrected, lai)
  )

# Prepare light attenuation data ----
tower_files <- list.files("data_working/neon_flux", pattern="*.csv",
                          full.names=TRUE)

measure_heights <- amf_var_info()

# For each PAR observation we need to know the following:
# - What height was it measured at?
# - What was the PPFD value for this observation?
# - What was the TOC shortwave flux?
# - What proportion of incoming shortwave was diffuse?
cat("Preparing PAR attenuation data...")
tower_toc_rad <- lapply(tower_files, function(file) {
  # Extract site information from filename
  site_amf <- str_remove(basename(file), ".csv")
  site_neon <- site_info$site_neon[match(site_amf, site_info$site_ameriflux)]
  cat("\t", site_neon, "\n")
  
  # Open tower observations
  this_tower <- read_csv(file, col_types=cols())
  
  # Get measurement heights for this tower
  this_heights <- measure_heights %>% filter(Site_ID == site_amf)
  sw_vars <- grep("SW_IN_1_1_\\d", names(this_tower), value=TRUE)
  ppfd_vars <- grep("PPFD_IN_\\d_\\d_\\d", names(this_tower), value=TRUE)
  ppfd_z <- this_heights %>%
    filter(Variable %in% ppfd_vars) %>% 
    select(Variable, Height)
  
  
  this_tower %>%
    mutate(
      # Attach site information for later
      SITE_NEON = site_neon,
      # Compute top of canopy values
      SW_IN_TOC = rowMeans(this_tower %>% select(all_of(sw_vars)), 
                           na.rm=TRUE),
      PPFD_IN_TOC = PPFD_IN_1_1_1,
      proportion_diffuse = SW_DIF / SW_IN_TOC,
      proportion_diffuse = pmin(pmax(proportion_diffuse, 0), 1)
    ) %>%
    # Drop observations without TOC values and values at night
    filter(!is.nan(SW_IN_TOC), !is.na(PPFD_IN_TOC), SW_IN_TOC > 100) %>%
    # Drop unnecessary columns
    select(SITE_NEON, PPFD_IN_TOC, proportion_diffuse, 
           all_of(ppfd_vars)) %>%
    # Pivot PAR observations to join in measurement heights
    pivot_longer(all_of(ppfd_vars), 
                 names_to = "ppfd_variable", values_to="ppfd_value") %>%
    left_join(ppfd_z, by=c("ppfd_variable"="Variable")) %>%
    # Cleanup
    rename(ppfd_z = Height) %>%
    filter(!is.na(ppfd_value))
  
}) %>% bind_rows()

# Fit LiDAR attenuation constants ----
fit_lidar_constants <- function(site, site_lai, zmax=NA, zmin=1) {
  # Z where cumulative LAI == total LAI and measurement height of Io
  if (is.na(zmax)) {
    zmax <- tower_toc_rad %>% filter(SITE_NEON == site) %>% 
      pull(ppfd_z) %>% max() %>% ceiling()
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
  this_tower_toc_rad <- tower_toc_rad %>% filter(SITE_NEON == site)
  
  iv_io_df <- this_tower_toc_rad %>%
    mutate(iv = ppfd_value,
           io = PPFD_IN_TOC,
           proportion_diffuse_bin = round(proportion_diffuse, digits=1)) %>%
    mutate(iv_io = iv / io) %>%
    group_by(ppfd_z, proportion_diffuse_bin) %>%
    summarize(
      iv_io_mean = suppressWarnings(log(iv_io)) %>% median(na.rm=TRUE) %>% exp(),
      # https://en.wikipedia.org/wiki/Geometric_standard_deviation
      iv_io_sd = suppressWarnings(log(iv_io / iv_io_mean)^2) %>% 
        mean(na.rm=TRUE) %>% sqrt() %>% exp()
    ) %>%
    drop_na() %>%
    mutate(cum_lai = approx(lad_profile$z, lad_profile$cum_lai, ppfd_z, rule=2)$y,
           site = site)
  
  # Force the y-intercept at log(iv_io) = 0
  iv_io_lm <- lm(
    log(iv_io_mean) ~ cum_lai + cum_lai:proportion_diffuse_bin + 0, 
    data=iv_io_df
  )
  
  return(list(
    site=site,
    lad_profile=lad_profile,
    iv_io_data=iv_io_df,
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

# Extract results ----

# LAD profiles
lad_data <- lapply(lidar_constants, function(l) {
  l$lad_profile %>% mutate(site = l$site)
}) %>% bind_rows()

# Light attenuation data
iv_io_data <- lapply(lidar_constants, function(l) l$iv_io_data) %>% bind_rows()

# Attenuation constants
lidar_constants_df <- lapply(lidar_constants, function(l) {
  list(
    site=l$site,
    lai_coef=l$lai_coef,
    # Force attenuation coefficients less than 0
    kb_real=min(unname(l$rad_coef[1]), 0),
    kd_real=min(unname(l$rad_coef[1] + l$rad_coef[2]), 0)
  )
}) %>% bind_rows()

# Plotting ----
use_sites <- c("ABBY", "DEJU", "JERC", "OSBS", "RMNP", "TALL", "WREF")
# Extrapolate lines for beam/diffuse coefficients
clai <- seq(0, 7, length.out=20)
iv_io_extrapolate <- expand.grid(clai, 1:nrow(lidar_constants_df)) %>%
  cbind(lidar_constants_df[.$Var2, ]) %>%
  select(-Var2) %>% rename(clai = Var1) %>%
  mutate(Beam = exp(kb_real * clai),
         Diffuse = exp(kd_real * clai)) %>%
  pivot_longer(c(Beam, Diffuse))

# Plot mean transmission points over the fitted line
(plot_mean_transmission <- iv_io_data %>%
  filter(site %in% use_sites) %>%
  ggplot(aes(x=cum_lai, y=iv_io_mean)) +
  geom_point(aes(color=proportion_diffuse_bin)) +
  geom_line(data=iv_io_extrapolate %>% filter(site %in% use_sites), 
            mapping=aes(x=clai, y=value, group=name, linetype=name)) +
  facet_wrap(~ site, scales="free") +
  labs(x="Cumulative LAI", y="Proportion available light",
       color="Proportion diffuse light", linetype="") +
  scale_y_log10() + theme_bw())

write_if_not_exist(lidar_constants_df, "data_out/neon_lidar_constants.csv")
write_if_not_exist(iv_io_data, "data_out/neon_lidar_transmission_data.csv")

# Plot LAD profiles for each site
(plot_lad_profiles <- ggplot(lad_data, aes(x=z)) +
  geom_line(aes(y=lad)) +
  coord_flip() +
  facet_wrap(~ site, scales="free") +
  labs(y=expression("Leaf area density"~frac(m^2, m^3)),
       x="Height above ground (m)") +
  theme_bw())

(plot_cum_lai_profiles <- ggplot(lad_data, aes(x=z)) +
    geom_line(aes(y=cum_lai)) +
    coord_flip() +
    facet_wrap(~ site, scales="free") +
    labs(y=expression("Cumulative LAI"~frac(m^2, m^2)),
         x="Height above ground (m)") +
    theme_bw())

lad_data %>%
  drop_na() %>%
  write_if_not_exist("data_out/neon_lad_profiles.csv")
