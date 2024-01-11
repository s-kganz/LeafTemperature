# To disaggregate whole-canopy photosynthesis, we assume that each slice's
# contribution to the whole-canopy value is proportional to the radiation
# flux reaching the slice. But, PS saturates beyond a certain point. Estimate
# this threshold at each site by comparing TOC flux with whole-canopy A.

library(tidyverse)
library(photosynthesis)

# Read/prepare data ----
flux <- read_csv("data_out/cross_site_flux_qc.csv") %>%
  select(site, timeBgn, data.fluxCo2.nsae.flux) %>%
  # Discard unreasonable values
  filter(data.fluxCo2.nsae.flux > -50, data.fluxCo2.nsae.flux < 20)

tower <- read_csv("data_out/cross_site_tc_data.csv") %>%
  filter(z == zmax) %>%
  select(TIMESTAMP, SITE_NEON, SW_IN_TOC, PPFD_IN_TOC, USTAR)

flux_tower_join <- flux %>%
  inner_join(tower,
             by=c("timeBgn"="TIMESTAMP", "site"="SITE_NEON")) %>%
  # Very strange values when shortwave flux is near zero, drop.
  # Observations when friction velocity is low can be biased, drop.
  filter(USTAR > 0.4, SW_IN_TOC > 50) %>%
  # Model fitting only works when the curve is downward-opening. Invert the
  # flux value to account for that.
  mutate(A = -data.fluxCo2.nsae.flux)

# Plot raw data ----
ggplot(flux_tower_join, aes(x=SW_IN_TOC, y=A)) +
  geom_point(alpha=0.1) +
  facet_wrap(~ site, scales="free_y") +
  labs(x="Top-of-canopy shortwave (W / m2)",
       y="Canopy photosynthesis (umol CO2 / m2 / s")

# Fit models ----
fits <- flux_tower_join %>%
  split(~ site) %>%
  map(fit_photosynthesis, .photo_fun="aq_response",
      .vars=list(.A=A, .Q=SW_IN_TOC), .method="ls") %>%
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
  geom_point(mapping=aes(x=SW_IN_TOC, y=A), data=flux_tower_join, alpha=0.1) +
  geom_line(mapping=aes(x=Q, y=A, color="Fit"), data=fit_pred, lwd=1) +
  facet_wrap(~ site, scales="free_y") +
  labs(x="Top-of-canopy shortwave (W / m2)",
       y="Canopy photosynthesis (umol CO2 / m2 / s)",
       color="") +
  theme_bw()

# Save model coefficients ----
write_csv(fits, "data_out/cross_site_aq_constants.csv")
