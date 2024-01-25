# This file calculates canopy-scale g1 for the Medlyn stomatal conductance
# model.

library(tidyverse)
library(bigleaf)
library(RcppRoll)
source("scripts/energy_balance/air.R")
source("scripts/energy_balance/leaf.R")

fit_medlyn_model_parameters <- function(site, tower, flux, zr, zh, L) {
  cat("Now processing", site, "\n")
  # Select needed columns from tower ----
  
  # For CO2, the MIXING_RATIO name means the mole fraction is per mole of dry 
  # air, while the CO2 name means the mole fraction is per mole of wet air. 
  # Both are in ppm and I don't think the differences matters so just lump 
  # the two versions together.
  co2_regex <- "CO2.+1_1_\\d"
  sw_regex  <- "SW_IN_1_1_\\d"
  ta_regex  <- "TA_1_1_\\d"
  g_regex   <- "G_1_1_\\d"
  
  # Other vars that don't need to be averaged. Don't pull LE because we use the
  # quality controlled fluxes from NEON.
  req_vars <-
    c("TIMESTAMP",
      "WS_1_1_1",
      "USTAR",
      "P",
      "PA",
      "RH",
      "H",
      "LW_IN",
      "LW_OUT",
      "SW_OUT")
  
  tower_select <- tower %>%
    select(all_of(req_vars)) %>%
    mutate(
      SW_IN = rowMeans(tower %>% select(matches(sw_regex)), na.rm = TRUE),
      TA = rowMeans(tower %>% select(matches(ta_regex)), na.rm = TRUE),
      CO2 = rowMeans(tower %>% select(matches(co2_regex)), na.rm = TRUE),
      G = rowMeans(tower %>% select(matches(g_regex)), na.rm = TRUE),
      NETRAD = SW_IN + LW_IN - SW_OUT - LW_OUT,
      VPD = vapor_pressure_deficit(TA + 273.15, RH)
    ) %>%
    rename(WS = WS_1_1_1) %>% drop_na()
  
  cat("Initial # observations:", nrow(tower_select), "\n")
  
  # Identify periods where T/ET is close to 1 ----
  no_precip_72hr <- tower_select %>%
    mutate(date = as.Date(TIMESTAMP)) %>%
    group_by(date) %>%
    summarize(P = sum(P, na.rm = TRUE)) %>%
    mutate(P_roll_sum = roll_sumr(P, n = 4)) %>%
    filter(P_roll_sum == 0)
  
  tower_transp_only <- tower_select %>%
    filter(as.Date(TIMESTAMP) %in% no_precip_72hr$date) %>%
    # Drop night time data
    filter(SW_IN > 100)
  
  cat("...after T/ET filtering", nrow(tower_transp_only), "\n")
  
  # Calculate roughness parameters ----
  roughness <- roughness.parameters(
    method = "wind_profile",
    zh = zh,
    zr = zr,
    Tair = tower_transp_only$TA,
    pressure = tower_transp_only$PA,
    wind = tower_transp_only$WS,
    ustar = tower_transp_only$USTAR,
    H = tower_transp_only$H
  )
  
  # Join QC'd flux measurements ----
  tower_transp_only_flux <- tower_transp_only %>%
    inner_join(flux, by = c("TIMESTAMP" = "timeBgn")) %>%
    rename(LE = data.fluxH2o.nsae.flux)
  
  cat("...after joining flux", nrow(tower_transp_only_flux), "\n")
  
  if (nrow(tower_transp_only_flux) < 50) {
    cat("Too few observations! Returning NA...\n")
    return(list(site=site, g0=NA, g1=NA, data=NA))
  }
  
  # Calculate conductances ----
  tower_gs <- tower_transp_only_flux %>%
    bind_cols(
      # Calculate Ga_h and Gb_h
      aerodynamic.conductance(
        .,
        .$TA,
        .$PA,
        .$WS,
        .$USTAR,
        .$H,
        zr = zr,
        zh = zh,
        d = roughness$d,
        z0m = roughness$z0m,
        # Dl is leaf characteristic dimension - we are working with needleleaves
        # here and other papers use this estimate
        Dl = 0.01,
        wind_profile = TRUE
      )
    ) %>%
    mutate(
      # Calculate corrected air T, CO2 conc, VPD after Knauer et al. (2018). These
      # seem to be less accurate than top-of-canopy measurements based on
      # instrumentation lower down on the tower so I don't use them.
      rho = dry_air_density(TA + 273.15, PA),
      gamma = psychrometric_constant(PA),
      e = saturation_vapor_pressure(TA + 273.15) * (1 - RH / 100),
      TA_corr = TA + H / (rho * Ga_h * eb_constants_$cp),
      e_corr = e + (LE * gamma) / (rho * Ga_h * eb_constants_$cp),
      VPD_corr = saturation_vapor_pressure(TA_corr + 273.15) - e_corr,
      CO2_corr = CO2 + data.fluxCo2.nsae.flux / Ga_h
    ) %>%
    bind_cols(
      # Finally, calculate surface conductance!
      surface.conductance(
        .,
        Tair = .$TA,
        pressure = .$PA,
        Rn = .$NETRAD,
        G = .$G,
        S = 0,
        VPD = .$VPD,
        LE = .$LE,
        Ga = .$Ga_h
      )
    ) %>%
    filter(
      # Discard NAs, negatives
      Gs_mol >= 0, !is.na(Gs_mol)
    )
    
  # Fit Medlyn model ----

  # A couple options for fitting.
  
  # V1: numerical optimization. This is more sensitive to outliers, mitigate by
  # weighting by sqrt(VPD)?
  medlyn_optimizer_helper <- function(Gs, VPD, GPP, CO2) {
    optimize(function(g1) {
      sum(((medlyn_gs(GPP, VPD, CO2, g1) - Gs) ^ 2) * sqrt(VPD))
    },
    c(1, 10))
  }
  
  fit_v1 <-
    with(tower_gs, medlyn_optimizer_helper(Gs_mol, VPD, GPP_U50_f, CO2))
  
  # V2: Theil-Sen. This option is less sensitive to outliers but it doesn't
  # strictly follow the model specification.
  # fit_v2 <- theil_sen_regression(Gs_mol ~ GPP / (sqrt(VPD) * CO2), data=tower_gs)
  
  tower_gs$Gs_mol_predicted <- medlyn_gs(
    tower_gs$GPP_U50_f, tower_gs$VPD, tower_gs$CO2, fit_v1$minimum
  )
  tower_gs$site <- site
  
  list(site=site,
       g0 = 0,
       g1 = fit_v1$minimum,
       data = tower_gs)
}

flux <- read_csv("data_out/cross_site_flux_partition_qc.csv")
lai  <- read_csv("data_out/neon_sampled_dhp_lai.csv")
meta <- read_csv("data_working/neon_site_metadata.csv") %>%
  left_join(lai, by=c("site_neon"="site")) %>%
  mutate(L_best = ifelse(is.na(lai), L_dhp, lai))

medlyn_fits <- lapply(1:nrow(meta), function(i) {
  this_site_amf  <- meta$site_ameriflux[i]
  this_site_neon <- meta$site_neon[i]
  this_zr        <- meta$tower_height[i]
  this_zh        <- meta$canopy_height[i]
  this_L         <- meta$L_best[i]
  
  this_flux <- flux %>%
    filter(site == this_site_neon)
  
  if (nrow(this_flux) == 0) {
    cat(this_site_neon, "has no flux! Skipping...\n")
    return(list(site=this_site_amf, g0=NA, g1=NA, data=NA))
  }
  
  this_tower_file <- file.path("data_working/neon_flux", str_c(this_site_amf, ".csv"))
  
  this_tower <- read_csv(this_tower_file, col_types=cols()) # Quiet, you!
  
  fit_medlyn_model_parameters(this_site_amf, this_tower,
                              this_flux, this_zr, this_zh, this_L)
})

medlyn_data <- lapply(medlyn_fits, function(x) x$data)
medlyn_df <- medlyn_data[!is.na(medlyn_data)] %>% bind_rows()
medlyn_coefs <- lapply(medlyn_fits, function(x) x[-4]) %>% bind_rows()

ggplot(NULL) + 
  geom_point(mapping=aes(x=Gs_mol, y=Gs_mol_predicted, color=VPD),
             data=medlyn_df) +
  geom_abline(slope=1, intercept=0, color="gray") +
  geom_label(mapping=aes(label=str_c("g1 = ", round(g1, digits=2))),
             x=0.4,
             y=1.5,
             data=medlyn_coefs) +
  facet_wrap(~ site) +
  xlim(0, 1.5) + coord_equal() +
  scale_color_viridis_c() +
  theme_bw() +
  labs(y = expression("Model G"[s] ~ bgroup("(", frac(mol, m ^ 2 ~ s), ")")),
       x = expression("Penman-Monteith G"[s] ~ bgroup("(", frac(mol, m ^
                                                                  2 ~ s), ")")),
       color = "VPD (kPa)")

write_csv(medlyn_df, "data_out/cross_site_medlyn_fit_results.csv")
write_csv(medlyn_coefs, "data_out/cross_site_medlyn_coefficients.csv")
