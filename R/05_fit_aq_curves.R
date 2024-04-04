# To disaggregate whole-canopy photosynthesis, we assume that each slice's
# contribution to the whole-canopy value is proportional to the radiation
# flux reaching the slice. But, PS saturates beyond a certain point. Estimate
# this threshold at each site by comparing TOC flux with whole-canopy A.

library(tidyverse)
library(photosynthesis)
#source("scripts/tower_util.R")

fit_aq_curves <- function(site_meta, site_flux_qc, tower_dir, outdir) {
  
  tower_files <- list.files(tower_dir, pattern="*.csv",
                            full.names=TRUE)
  
  tower_toc_rad <- lapply(tower_files, function(file) {
    site_amf <- str_remove(basename(file), ".csv")
    site_neon <- site_meta$site_neon[match(site_amf, site_meta$site_ameriflux)]
    this_tower <- read_csv(file, col_types=cols())
    
    this_tower %>%
      mutate(
        SITE_NEON = site_neon,
        SW_IN_TOC = rowMeans(this_tower %>% select(matches("SW_IN_1_1_\\d")), 
                             na.rm=TRUE)
      ) %>%
      select(TIMESTAMP, SITE_NEON, SW_IN_TOC)
  }) %>% bind_rows()
  
  flux_tower_join <- site_flux_qc %>%
    inner_join(tower_toc_rad,
               by=c("timeBgn"="TIMESTAMP", "site"="SITE_NEON")) %>%
    # Very strange values when shortwave flux is near zero, drop.
    # Observations when friction velocity is low can be biased, drop.
    filter(SW_IN_TOC > 50) %>%
    mutate(.A = -data.fluxCo2.nsae.flux,
           .Qabs = SW_IN_TOC)
  
  # Plot raw data ----
  # ggplot(flux_tower_join, aes(x=.Qabs, y=.A)) +
  #   geom_point(alpha=0.1) +
  #   facet_wrap(~ site, scales="free_y") +
  #   labs(x="Top-of-canopy shortwave (W / m2)",
  #        y="Canopy photosynthesis (umol CO2 / m2 / s") +
  #   theme_bw()
  
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
  
  # Save model coefficients and fit ----
  write_if_not_exist(
    fits, 
    file.path(outdir, "cross_site_aq_constants.csv")
  )
  write_if_not_exist(
    flux_tower_join, 
    file.path(outdir, "cross_site_aq_data.csv")
  )
}
