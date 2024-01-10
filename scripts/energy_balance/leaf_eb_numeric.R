# Another version of leaf energy balance that must be solved numerically.

source("scripts/energy_balance/leaf.R")
source("scripts/energy_balance/radiation.R")
source("scripts/energy_balance/air.R")

#' Numerical leaf energy balance model
#'
#' @param Tl Leaf temperature (K)
#' @param Ta Air temperature (K)
#' @param u Wind speed (m/s)
#' @param gs Stomatal conductance (mol / m2 sec)
#' @param Pa Air pressure (kPa) 
#' @param RH Relative humidity (%)
#' @param SW_abs Absorbed shortwave radiation (W/m2)
#' @param LW_abs_dn Absorbed downwelling longwave radiation (W/m2)
#' @param a_lw Longwave absorptance (0-1)
#' @param a_sw Shortwave absorptance (0-1)
#' @param d Leaf characteristic dimension (m)
#'
#'
#' @return
#' @export
#'
#' @examples
energy_balance_error <- function(Tl, Ta, u, gs, Pa, RH, 
                                 SW_abs, LW_abs_dn,
                                 a_lw=0.98, a_sw=0.50, d=0.01) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  vpd     <- vapor_pressure_deficit(Ta, RH)
  gamma   <- psychrometric_constant(Pa, eb_constants_$cp, eb_constants_$MWr, 
                                    eb_constants_$l)
  desat   <- esat_slope(Ta)
  
  # Conductances
  gbH   <- conductance_boundary_heat_needleleaf(u, d, rho_mol)
  gR    <- conductance_radiative_heat(Ta, Pa, rho, rho_mol, a_lw, 
                                      eb_constants_$cp)
  gtot  <- total_conductance(gs, gbH)
  
  # Net radiation
  Rn <- net_nonisothermal_radiation(SW_abs, LW_abs_dn, Tl, a_lw)
  
  # Latent heat flux
  omega  <- decoupling_coefficient(gR, gbH, gs, desat, gamma, a_lw)
  LE_imp <- latent_heat_imposed(gtot, vpd, Pa, eb_constants_$l)
  LE_eq  <- latent_heat_equilibrium(Rn, desat, gamma, gR, gbH)
  LE     <- omega * LE_eq + (1 - omega) * LE_imp
  
  # Sensible heat flux
  H <- sensible_heat(gbH, Tl, Ta)
  
  # Final energy closure calculation
  abs(Rn - H - LE)
}

#' @describeIn energy_balance_error User-friendly driver function to run
#' numerical optimization.
energy_balance_driver <- function(Ta, u, gs, Pa, RH, 
                                  SW_abs, LW_abs_dn,
                                  a_lw=0.98, a_sw=0.50, d=0.01) {
  optimize(
    energy_balance_error,
    c(Ta-10, Ta+10),
    Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn,
    a_lw=a_lw, a_sw=a_sw, d=d,
    tol=1 # solve to within 1 W of closure
  )
}

if (sys.nframe() == 0) {
  # Test optimization
  Ta <- 298
  u <- 2
  gs <- 0.01
  Pa <- 100
  RH <- 0.5
  SW_abs <- 400
  LW_abs_dn <- 200
  
  Tl_test <- seq(Ta-10, Ta+10, length.out=50)
  eb_error <- energy_balance_error(Tl_test, Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn)
  
  optim_result <- energy_balance_driver(
    Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn
  )
  
  plot(Tl_test, eb_error, col="black")
  abline(h=optim_result$objective, col="red", lty="dashed")
  abline(v=optim_result$minimum, col="red", lty="dashed")
}
