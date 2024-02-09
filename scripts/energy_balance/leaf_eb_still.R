source("scripts/energy_balance/air.R")
source("scripts/energy_balance/radiation.R")
source("scripts/energy_balance/leaf.R")
source("scripts/energy_balance/constants.R")

# This EB function is a modified version of the approach in Still et al. (2022).
# In this version we use leaf temperature to calculate leaf-air VPD and net
# radiation instead of substituting in air temperature. This means that we have
# to solve the EB equation numerically.
leaf_temperature_error_function <- function(
    Tl, Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn,
    a_lw=0.98, a_sw=0.50, d=0.01
) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  # VPD is calculated from the leaf to the air
  e_leaf  <- saturation_vapor_pressure(Tl)
  e_air   <- saturation_vapor_pressure(Ta) * RH/100
  
  vpd     <- max(e_leaf - e_air, 0)
  gamma   <- psychrometric_constant(Pa, eb_constants_$cp, eb_constants_$MWr, 
                                    eb_constants_$l)
  desat   <- esat_slope(Ta)
  
  # Conductances
  gbH   <- conductance_boundary_heat_needleleaf(u, d, rho_mol)
  gR    <- conductance_radiative_heat(Ta, Pa, rho, rho_mol, 
                                      a_lw, eb_constants_$cp)
  gtot  <- total_conductance(gs, gbH)
  
  # Net radiation - but using Tl instead of Ta
  Rn <- net_isothermal_radiation(Sw_abs, LW_abs_dn, Tl, a_lw)
  
  # Latent heat fluxes
  omega  <- decoupling_coefficient(gR, gbH, gs, desat, gamma, a_lw)
  LE_imp <- latent_heat_imposed(gtot, vpd, Pa, eb_constants_$l)
  LE_eq  <- latent_heat_equilibrium(Rn, desat, gamma, gR, gbH)
  LE     <- omega * LE_eq + (1 - omega) * LE_imp
  
  # Final leaf temperature calcluation
  Y  <- gbH / (gbH + gR)
  Tl_calc <- Ta + Y * (Rn - LE) / (eb_constants_$cp * eb_constants_$Mair * gbH)
  
  return(abs(Tl - Tl_calc))
}

# Helper func to optimize the above EB model
leaf_temperature_nonisothermal <- function(Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn,
                                           a_lw=0.98, a_sw=0.50, d=0.01,
                                           bounds=20, aslist=FALSE) {
  tryCatch(
    {
      optimize(
        function(Tl) {leaf_temperature_error_function(
          Tl, Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn, a_lw, a_sw, d
        )},
        c(Ta - bounds, Ta + bounds)
      )
    },
    error=function(e) {
      # print(e)
      return(NA)
    }
  )
}

# This EB function is equivalent to the procedure in Still et al. (2022). It is
# called isothermal because air temperature substitutes for leaf temperature in
# calculating net radiation and VPD. This does not use optimization to solve
# for leaf T.
leaf_temperature_isothermal <- function(Ta, u, gs, Pa, RH, SW_abs, LW_abs_dn,
                                        a_lw=0.98, a_sw=0.50, d=0.01,
                                        aslist=FALSE) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, eb_constants_$cp, eb_constants_$Mair)
  vpd     <- vapor_pressure_deficit(Ta, RH)
  gamma   <- psychrometric_constant(Pa, eb_constants_$cp, eb_constants_$MWr, 
                                    eb_constants_$l)
  desat   <- esat_slope(Ta)
  
  # Conductances
  gbH   <- conductance_boundary_heat_needleleaf(u, d, rho_mol)
  gR    <- conductance_radiative_heat(Ta, Pa, rho, rho_mol, 
                                      a_lw, eb_constants_$cp)
  gtot  <- total_conductance(gs, gbH)
  
  # Net isothermal radiation
  Rn <- net_isothermal_radiation(SW_abs, LW_abs_dn, Ta, a_lw)
  
  # Latent heat fluxes
  omega  <- decoupling_coefficient(gR, gbH, gs, desat, gamma, a_lw)
  LE_imp <- latent_heat_imposed(gtot, vpd, Pa, eb_constants_$l)
  LE_eq  <- latent_heat_equilibrium(Rn, desat, gamma, gR, gbH)
  LE     <- omega * LE_eq + (1 - omega) * LE_imp
  
  # Final leaf temperature calcluation
  Y  <- gbH / (gbH + gR)
  Tl <- Ta + Y * (Rn - LE) / (eb_constants_$cp * eb_constants_$Mair * gbH)
  
  # Temperature forcings calculation
  dT_Rn <- Y * Rn / (eb_constants_$cp * eb_constants_$Mair * gbH)
  dT_LE <- Y * LE / (eb_constants_$cp * eb_constants_$Mair * gbH)
  
  Tl <- Ta + dT_Rn - dT_LE
  
  if (aslist) {
    list(
      Tl=Tl,
      dT_Rn=dT_Rn,
      dT_LE=dT_LE,
      omega=omega,
      gtot=gtot,
      gbH=gbH,
      LE=LE,
      LE_imp=LE_imp,
      LE_eq=LE_eq,
      Rn=Rn
    )
  } else {
    return(Tl)
  }
}

if (sys.nframe() == 0) {
  # Example parameters
  leaf_temperature_isothermal(
    300, 2, 0.2, 100, 50, 100, 400
  )
}