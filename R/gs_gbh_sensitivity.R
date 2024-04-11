# This file looks at the relationship between gs, gbh, and the leaf-air 
# temperature difference.


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
  optim_result <- stats::optimize(
    energy_balance_conductance_error,
    c(Ta-bounds, Ta+bounds),
    Ta, gs, gbH, Pa, RH, SW_dn, LW_dn,
    a_lw=a_lw, a_sw=a_sw,
    tol=0.01 # solve to within 0.01 K of optimum
  )
  
  return(optim_result$minimum)
}

#' Stomatal/boundary layer sensitivity analysis
#'
#' @param outdir Directory where model results will be written.
#' @param gs_min Minimum value for stomatal conductance (mol m-2 s-1).
#' @param gs_max Maximum value for stomatal conductance (mol m-2 s-1).
#' @param gbH_min Minimum value for boundary layer conductance (mol m-2 s-1).
#' @param gbH_max Maximum value for boundary layer conductance (mol m-2 s-1).
#' @param Ta Air temperature (K).
#' @param Pa Air pressure (kPa).
#' @param RH Relative humidity (%, from 0-100).
#' @param a_lw Longwave absorptivity of leaves.
#' @param a_sw Shortwave absorptivity of leaves.
#'
#' @description
#' Stomatal conductance and boundary layer conductance are two important parameters that drive the temperature of a leaf. This function evaluates how these parameter drive temperature under a low, medium, and high radiation load.
#'
#' @export
#'
gs_gbh_sensitivity <- function(outdir, 
                               gs_min=0.01, gs_max=0.5, gbH_min=0.75, gbH_max=5,
                               Ta=293, Pa=100, RH=50, a_lw=0.98, a_sw=0.50) {
  # Set parameters ----
  gs_values <- seq(gs_min, gs_max, length.out=100)
  gbH_values <- seq(gbH_min, gbH_max, length.out=100)
  # Generally shortwave is about 1/3 of the absorbed radiation load
  # in the canopy top.
  tot_rad <- c(800, 1000, 1200)
  
  grid <- expand.grid(gs_values, gbH_values, tot_rad) %>%
    rename(
      gs = Var1,
      gbH = Var2,
      tot_rad = Var3
    ) %>%
    # Set defaults
    mutate(
      Ta = 293, # K
      Pa = 100, # kPa
      RH = 50, # %
      SW_dn = tot_rad * (1/3),
      LW_dn = tot_rad * (2/3),
      a_lw = 0.98, # -
      a_sw = 0.50  # -
    )
  
  # Run energy balance ----
  message("Running gs/gbh model...")
  grid_eb <- grid %>%
    mutate(
      Tl = purrr::pmap_dbl(list(Ta, gs, gbH, Pa, RH, SW_dn, LW_dn),
                           energy_balance_conductance_driver),
      dT = Tl - Ta
    )
  
  write_if_not_exist(grid_eb, file.path(outdir, "gs_gbh_sensitivity.csv"))
}
       