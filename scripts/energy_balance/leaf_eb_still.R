source("scripts/energy_balance/air.R")
source("scripts/energy_balance/radiation.R")
source("scripts/energy_balance/leaf.R")

#' Back-calculate stomatal conductance from leaf temperature
#' and surface conditions.
#'
#' @param envpar 
#' @param constants 
#'
#' @return
#' @export
#'
#' @examples
stomatal_conductance <- function(Tl, envpar, constants) {
  optim_result <- optimize(
    function(gs) {
      envpar[["gs"]] <- gs
      abs(Tl - leaf_temperature(envpar, constants)$Tl)
    },
    c(0, 1) # optimization limits
  )
  
  optim_result$minimum
}

# This EB function is a modified version of the approach in Still et al. (2022).
# In this version we use leaf temperature to calculate leaf-air VPD and net
# radiation instead of substituting in air temperature. This means that we have
# to solve the EB equation numerically.
leaf_temperature_error_function <- function(
    Tl, Ta, Pa, RH, u, gs, SW_dir, SW_dif, SW_out, LW_down, G, constants
) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, constants$cp, constants$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, constants$cp, constants$Mair)
  # VPD is calculated from the leaf to the air
  e_leaf  <- saturation_vapor_pressure(Tl)
  e_air   <- saturation_vapor_pressure(Ta) * RH/100
  
  vpd     <- max(e_leaf - e_air, 0)
  gamma   <- psychrometric_constant(Pa, constants$cp, constants$MWr, constants$l)
  desat   <- esat_slope(Ta)
  
  # Conductances
  gbH   <- conductance_boundary_heat_needleleaf(u, constants$d, rho_mol)
  gR    <- conductance_radiative_heat(Ta, Pa, rho, rho_mol, 
                                      constants$a_lw, constants$cp)
  gtot  <- total_conductance(gs, gbH)
  
  # Net radiation - but using Tl instead of Ta
  Rn <- net_isothermal_radiation(SW_dir, SW_dif, SW_out, LW_down, Tl, G, 
                                 constants$a_sw, constants$a_lw)
  
  # Latent heat fluxes
  omega  <- decoupling_coefficient(gR, gbH, gs, desat, gamma, constants$a_lw)
  LE_imp <- latent_heat_imposed(gtot, vpd, Pa, constants$l)
  LE_eq  <- latent_heat_equilibrium(Rn, desat, gamma, gR, gbH)
  LE     <- omega * LE_eq + (1 - omega) * LE_imp
  
  # Final leaf temperature calcluation
  Y  <- gbH / (gbH + gR)
  Tl_calc <- Ta + Y * (Rn - LE) / (constants$cp * constants$Mair * gbH)
  
  return(abs(Tl - Tl_calc))
}

# Helper func to optimize the above EB model
leaf_temperature_nonisothermal <- function(Ta, Pa, RH, u, gs, SW_dir, SW_dif,
                                           SW_out, LW_down, G, constants) {
  tryCatch(
    {
      optimize(
        function(Tl) {leaf_temperature_error_function(
          Tl, Ta, Pa, RH, u, gs, SW_dir, SW_dif, SW_out, LW_down, G, constants
        )},
        c(Ta - 10, Ta + 10)
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
leaf_temperature_isothermal <- function(Ta, Pa, RH, u, gs, SW_dir, SW_dif, 
                                        SW_out, LW_down, G, constants) {
  # Physical parameters of the air
  rho     <- dry_air_density(Ta, Pa, constants$cp, constants$Mair)
  rho_mol <- dry_air_molar_density(Ta, Pa, constants$cp, constants$Mair)
  vpd     <- vapor_pressure_deficit(Ta, RH)
  gamma   <- psychrometric_constant(Pa, constants$cp, constants$MWr, constants$l)
  desat   <- esat_slope(Ta)
  
  # Conductances
  gbH   <- conductance_boundary_heat_needleleaf(u, constants$d, rho_mol)
  gR    <- conductance_radiative_heat(Ta, Pa, rho, rho_mol, 
                                      constants$a_lw, constants$cp)
  gtot  <- total_conductance(gs, gbH)
  
  # Net isothermal radiation
  Rn <- net_isothermal_radiation(SW_dir, SW_dif, SW_out, LW_down, Ta, G, 
                      constants$a_sw, constants$a_lw)
  
  # Latent heat fluxes
  omega  <- decoupling_coefficient(gR, gbH, gs, desat, gamma, constants$a_lw)
  LE_imp <- latent_heat_imposed(gtot, vpd, Pa, constants$l)
  LE_eq  <- latent_heat_equilibrium(Rn, desat, gamma, gR, gbH)
  LE     <- omega * LE_eq + (1 - omega) * LE_imp
  
  # Final leaf temperature calcluation
  Y  <- gbH / (gbH + gR)
  Tl <- Ta + Y * (Rn - LE) / (constants$cp * constants$Mair * gbH)
  
  list(
    Tl=Tl,
    omega=omega,
    gtot=gtot,
    gbH=gbH,
    LE=LE,
    LE_imp=LE_imp,
    LE_eq=LE_eq,
    Y=Y,
    Rn=Rn
  )
}

if (sys.nframe() == 0) {
  # Example parameters
  eb_envpar <- list(
    Ta=298,        # Air temperature, K
    Pa=98,         # Air pressure, kPa
    RH=60,         # Relative humidity, % (0 - 100)
    u=2,           # Wind speed m/sec
    gs=0.4,        # Stomatal conductance mol/m2 sec
    LW_down=350,   # Downwelling longwave radiation W/m2
    SW_dir=700,    # Downwelling direct shortwave radiation W/m2
    SW_dif=150,    # Downwelling diffuse shortwave ardiation W/m2
    SW_out=60,     # Upwelling shortwave radiation
    is_sunlit=1,   # Is this leaf in shade? (1/0)
    illuminance=1  # Relative illuminance of the leaf (0-1)
  )
}