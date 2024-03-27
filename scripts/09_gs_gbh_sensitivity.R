# This file looks at the relationship between gs, gbh, and the leaf-air 
# temperature difference.

library(tidyverse)
library(colorspace)
library(metR)
source("scripts/energy_balance/leaf_eb_numeric.R")
source("scripts/tower_util.R")


# Driver functions ----
# Have to write separate EB function to manipulate conductances directly.
# TODO consider refactoring original version to make sensitivity analysis
# like this easier.
energy_balance_conductance_error <- function(Tl, Ta, gs, gbH, Pa, RH,
                                             SW_dn, LW_dn,
                                             a_lw=0.98, a_sw=0.50) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  vpd     <- vapor_pressure_deficit(Tl, RH)
  gamma   <- psychrometric_constant(Pa, eb_constants_$cp, eb_constants_$MWr, 
                                    eb_constants_$l)
  desat   <- esat_slope(Ta)
  
  # Conductances - don't need to calculate gbH since that is provided
  gR    <- conductance_radiative_heat(Ta, Pa, rho, rho_mol, a_lw, 
                                      eb_constants_$cp)
  gtot  <- total_conductance(gs, gbH)
  
  # Net radiation - here we don't consider canopy effects, so assume
  # L = 0.
  SW_abs <- SW_dn * a_sw
  LW_abs <- LW_dn * a_lw
  Rn <- net_nonisothermal_radiation(SW_abs, LW_abs, Tl, a_lw)
  
  # Latent heat flux
  omega  <- decoupling_coefficient(gR, gbH, gs, desat, gamma, a_lw)
  LE_imp <- latent_heat_imposed(gtot, vpd, Pa, eb_constants_$l)
  LE_eq  <- latent_heat_equilibrium(Rn, desat, gamma, gR, gbH)
  LE     <- omega * LE_eq + (1 - omega) * LE_imp
  
  # Sensible heat flux
  H <- sensible_heat(gbH, Tl, Ta)
  
  # Final energy closure calculation
  error <- abs(Rn - H - LE)
}

energy_balance_conductance_driver <- function(Ta, gs, gbH, Pa, RH, SW_dn, LW_dn,
                                              a_lw=0.98, a_sw=0.50, bounds=20) {
  optim_result <- optimize(
    energy_balance_conductance_error,
    c(Ta-bounds, Ta+bounds),
    Ta, gs, gbH, Pa, RH, SW_dn, LW_dn,
    a_lw=a_lw, a_sw=a_sw,
    tol=0.01 # solve to within 0.01 K of optimum
  )
  
  return(optim_result$minimum)
}

# Set parameters ----
gs_min <- 0.01
gs_max <- 0.5
gbH_min <- 0.75
gbH_max <- 5
gs_values <- seq(gs_min, gs_max, length.out=100)
gbH_values <- seq(gbH_min, gbH_max, length.out=100)

grid <- expand.grid(gs_values, gbH_values) %>%
  rename(
    gs = Var1,
    gbH = Var2
  ) %>%
  # Set defaults
  mutate(
    Ta = 293, # K
    Pa = 100, # kPa
    RH = 50, # %
    SW_dn = 800, # W m-2
    LW_dn = 600, # W m-2
    a_lw = 0.98, # -
    a_sw = 0.50  # -
  )

# Run energy balance ----
grid_eb <- grid %>%
  mutate(
    Tl = pmap_dbl(list(Ta, gs, gbH, Pa, RH, SW_dn, LW_dn),
                  energy_balance_conductance_driver,
                  .progress="EB conductance sensitivity"),
    dT = Tl - Ta
  )

write_if_not_exist(grid_eb, "data_out/model_runs/gs_gbh_sensitivity.csv")

# Plotting ----
# example gbH values
rho_mol <- dry_air_molar_density(293, 100)
eg_gbH <- data.frame(
  label=c("NL, windy", "NL, still", "BL, windy", "BL, still"),
  gbH=c(0.0060 * rho_mol * (5)^0.6 * (0.01)^-0.4,
        0.0060 * rho_mol * (0.5)^0.6 * (0.01)^-0.4,
        0.0105 * rho_mol * sqrt(5 / 0.1),
        0.0105 * rho_mol * sqrt(0.5 / 0.1))
)

# NEON model output
neon_eb <- read_csv("data_out/model_runs/cross_site_eb.csv")

base_contour_plot <- grid_eb %>%
  ggplot(aes(x=gs, y=gbH)) +
  geom_contour_filled(aes(z=dT)) +
  scale_fill_discrete_divergingx() +
  theme_bw() +
  labs(x=expression("g"[s]~"(mol m"^-2~"s"^-1*")"),
       y=expression("g"[bH]~"(mol m"^-2~"s"^-1*")"),
       fill=expression(Delta*"T (K)"))

base_contour_plot +
  geom_hline(data=eg_gbH, mapping=aes(yintercept=gbH), color="grey20",
             linetype="dashed") +
  geom_text(data=eg_gbH, mapping=aes(y=gbH, label=label), x=0.02,
            hjust=0, vjust=0, fontface="italic", color="grey20")

neon_eb_filter <- neon_eb %>%
  # Discard extremes
  filter(PS_LAYER_GS > gs_min, PS_LAYER_GS < gs_max,
         EB_MODEL_gbH > gbH_min, EB_MODEL_gbH < gbH_max)

ggplot(NULL) +
  geom_point(data=neon_eb_filter, mapping=aes(x=PS_LAYER_GS, y=EB_MODEL_gbH),
             alpha=0.1) +
  geom_contour(data=grid_eb, 
               mapping=aes(x=gs, y=gbH, z=dT),
               lwd=1.7, color="black") +
  geom_contour(data=grid_eb, 
               mapping=aes(x=gs, y=gbH, z=dT, color=after_stat(level)),
               lwd=1.5) +
  scale_color_divergent() +
  geom_label_contour(data=grid_eb, mapping=aes(x=gs, y=gbH, z=dT),
                     skip=0, label.placer=label_placer_fraction()) +
  theme_bw() +
  labs(x=expression("Stomatal conductance"~"(mol m"^-2~"s"^-1*")"),
       y=expression("Boundary layer heat conductance"~"(mol m"^-2~"s"^-1*")"),
       color=expression(Delta*"T (K)"))

neon_eb %>%
  mutate(EB_MODEL_dT = EB_MODEL_Tl - LAYER_TA - 273.15) %>%
  ggplot(aes(x=EB_MODEL_LE, y=EB_MODEL_H, color=EB_MODEL_dT)) +
  geom_point(alpha=0.5) +
  scale_color_divergent(mid="#f5f0d5") +
  geom_abline(slope=1, intercept=0, color="red") +
  theme_bw() +
  labs(x=expression("Latent heat flux (W m"^-2*")"),
       y=expression("Sensible heat flux (W m"^-2*")"),
       color=expression(Delta*"T (K)"))

       