# To disaggregate whole-canopy photosynthesis, we assume that each slice's
# contribution to the whole-canopy value is proportional to the radiation
# flux reaching the slice. But, PS saturates beyond a certain point. Estimate
# this threshold at each site by comparing TOC flux with whole-canopy A.

library(tidyverse)
library(photosynthesis)
source("scripts/tower_util.R")

# Read/prepare data ----
flux <- read_csv("data_out/cross_site_flux_partition_qc.csv")
meta <- read_csv("data_working/neon_site_metadata.csv")

tower_files <- list.files("data_working/neon_flux", pattern="*.csv",
                          full.names=TRUE)

tower_toc_rad <- lapply(tower_files, function(file) {
  site_amf <- str_remove(basename(file), ".csv")
  site_neon <- meta$site_neon[match(site_amf, meta$site_ameriflux)]
  this_tower <- read_csv(file, col_types=cols())
  
  this_tower %>%
    mutate(
      SITE_NEON = site_neon,
      SW_IN_TOC = rowMeans(this_tower %>% select(matches("SW_IN_1_1_\\d")), 
                           na.rm=TRUE)
    ) %>%
    select(TIMESTAMP, SITE_NEON, SW_IN_TOC)
}) %>% bind_rows()

flux_tower_join <- flux %>%
  inner_join(tower_toc_rad,
             by=c("timeBgn"="TIMESTAMP", "site"="SITE_NEON")) %>%
  # Very strange values when shortwave flux is near zero, drop.
  # Observations when friction velocity is low can be biased, drop.
  filter(SW_IN_TOC > 50) %>%
  mutate(.A = -data.fluxCo2.nsae.flux,
         .Qabs = SW_IN_TOC)

# Plot raw data ----
ggplot(flux_tower_join, aes(x=.Qabs, y=.A)) +
  geom_point(alpha=0.1) +
  facet_wrap(~ site, scales="free_y") +
  labs(x="Top-of-canopy shortwave (W / m2)",
       y="Canopy photosynthesis (umol CO2 / m2 / s") +
  theme_bw()

# Fit models ----
fits <- flux_tower_join %>%
  split(~ site) %>%
  map(
    # The default parameter bounds lead to numerical instability so we fit the
    # model manually with better constraints. The defaults seem to allow for
    # cases where the data doesn't include the asymptote. We always have the
    # asymptote so we can constrain a little more.
    function(.data) {
      minpack.lm::nlsLM(
        data = .data, 
        .A ~ marshall_biscoe_1980(Q_abs = .data[[".Qabs"]], k_sat, phi_J, theta_J) - Rd,
        # Attempt to estimate starting parameters
        start = photosynthesis:::get_init_aq_response(.data),
        # Set lower limits
        lower = c(min(.data[[".A"]]), 0, 0, 0),
        # set upper limits
        upper = c(2 * max(.data[[".A"]]), 0.5, 0.99, 5),
        # set max iterations for curve fitting
        control = nls.lm.control(maxiter = 100)
      )
    }
  ) %>%
  map(coef) %>%
  map(t) %>%
  map(as.data.frame) %>%
  imap_dfr(~ mutate(.x, site = .y))

# Show fit on top of raw data ----
sw_range <- seq(0, 1250)
fit_pred <- fits %>%
  mutate(
    A_pred = pmap(
      list(k_sat, phi_J, theta_J, Rd),
      function(k_sat, phi_J, theta_J, Rd) {
        data.frame(
          A = marshall_biscoe_1980(sw_range, k_sat, phi_J, theta_J) - Rd,
          Q = sw_range
        )
      }
    )
  ) %>%
  select(site, A_pred) %>% unnest(A_pred)

ggplot(NULL) +
  geom_point(mapping=aes(x=.Qabs, y=.A), data=flux_tower_join, alpha=0.1) +
  geom_line(mapping=aes(x=Q, y=A, color="Fit"), data=fit_pred, lwd=1) +
  facet_wrap(~ site, scales="free_y") +
  labs(x="Top-of-canopy shortwave (W / m2)",
       y="Canopy photosynthesis (umol CO2 / m2 / s)",
       color="") +
  theme_bw()

# Save model coefficients ----
write_if_not_exist(fits, "data_out/cross_site_aq_constants.csv")
