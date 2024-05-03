# Another version of leaf energy balance that must be solved numerically.

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
#' @param aslist Return error and fluxes as a list? If false, only returns error.
#' @param bounds Deviation from air temperature to define the search interval. For example, if Ta = 290 K and bounds = 10 K, the model will search for an optimal leaf temperature in the range of [280, 300] K.
#' @param run_no_transp Should the model be run a second time with gs = 0.001? This is necessary to calculate temperature forcings from the fluxes. If TRUE, additional outputs are returned.
#'
#' @return If aslist == TRUE, a list containing calculated heat fluxes, conductances, and decoupling. Otherwise, just the energy balance error.
#' @keywords internal
#'
#' 
energy_balance_error <- function(Tl, Ta, u, gs, Pa, RH, 
                                 SW_abs, LW_abs_dn,
                                 a_lw=0.95, a_sw=0.50, d=0.01,
                                 aslist=FALSE) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  vpd     <- vapor_pressure_deficit(Tl, RH)
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
  error <- abs(Rn - H - LE)
  
  # Return everything as a list
  if (aslist) {
    return(list(
      error=error,
      Rn=Rn,
      H=H,
      LE=LE,
      omega=omega,
      gbH=gbH,
      gR=gR,
      gtot=gtot
    ))
  } else {
    return(error)
  }
}

#' @describeIn energy_balance_error User-friendly driver function to run
#' numerical optimization.
energy_balance_driver <- function(Ta, u, gs, Pa, RH, 
                                  SW_abs, LW_abs_dn,
                                  a_lw=0.95, a_sw=0.50, d=0.01, bounds=20,
                                  run_no_transp=TRUE) {
  optim_result <- optimize(
    energy_balance_error,
    c(Ta-bounds, Ta+bounds),
    Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn,
    a_lw=a_lw, a_sw=a_sw, d=d,
    tol=0.01 # solve to within 0.01 K of optimum
  )
  
  Tl <- optim_result$minimum
  
  fluxes <- energy_balance_error(Tl, Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn,
                                 a_lw, a_sw, d, aslist=TRUE)
  
  if (run_no_transp) {
    optim_result_notransp <- optimize(
      energy_balance_error,
      c(Ta-bounds, Ta+bounds),
      Ta, u, 0, Pa, RH, SW_abs, LW_abs_dn,
      a_lw=a_lw, a_sw=a_sw, d=d,
      tol=0.01
    )
    
    Tl_notransp <- optim_result_notransp$minimum
    
    dT_Rn <- Tl_notransp - Ta
    dT_LE <- Tl_notransp - Tl
    
    return(c(Tl=Tl, fluxes, dT_Rn=dT_Rn, dT_LE=dT_LE))
    
  } else {
    
    return(c(Tl=Tl, fluxes))
  }
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
  abline(h=optim_result$error, col="red", lty="dashed")
  abline(v=optim_result$Tl, col="red", lty="dashed")
}
