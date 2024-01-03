#' Net isothermal radiation.
#'
#' @param SW_dir Direct downwelling shortwave radiation (W/m2)
#' @param SW_dif Diffuse downwelling shortwave radiation (W/m2)
#' @param SW_out Upwelling shortwave radiation (W/m2)
#' @param LW_down Downwelling longwave radiation (W/m2)
#' @param Ta Air temperature (K) 
#' @param is_sunlit Is the leaf in shade? (1/0)
#' @param illuminance Accessibility of diffuse light (0-1)
#' @param a_sw Absorptivity of shortwave radiation (0.5)
#' @param a_lw Absorptivity of longwave radiation (0.98)
#'
#' @return Rn, Net isothermal radiation (W/m2)
#' @export
#'
#' @examples
net_isothermal_radiation <- function(SW_dir, SW_dif, SW_out, LW_down, Ta, G,
                          a_sw=0.5, a_lw=0.98) {
  sb <- 5.67e-8
  SW_net <- SW_dir + SW_dif - SW_out
  Rn <- a_sw * SW_net + a_lw * (LW_down + G) - a_lw * sb * Ta^4
}

#' Density of dry air.
#' @param Ta Air temperature (K)
#' @param Pa Air pressure (kPa)
#' @param cp Heat capacity of dry air (J/kg K)
#' @param M Molar mass of dry air (kg / mol)
#'
#' @return rho, air density (kg / m3)
#' @export
#'
#' @examples
dry_air_density <- function(Ta, Pa, cp=1010, M=0.029) {
  R <- 8.314 # Univ gas constant
  Pa_pa <- Pa * 1000 # convert to Pa
  
  (Pa_pa * M) / (R * Ta)
}

#' Molar density of dry air.
#'
#' @param Ta Air temperature (K)
#' @param Pa Air pressure (Pa)
#' @param cp Heat capacity of dry air (J/kg K)
#' @param M Molar mass of dry air (kg/mol)
#'
#' @return rho_mol, molar density of dry air (mol / m3)
#' @export
#'
#' @examples
dry_air_molar_density <- function(Ta, Pa, cp=1010, M=0.029) {
  rho <- dry_air_density(Ta, Pa, cp, M) # kg/m3
  rho / M # kg / m3 * mol / kg = mol / m3
}

#' Leaf conductance to radiative heat transfer.
#'
#' @param Ta Air temperature (K)
#' @param Pa Air pressure (kPa)
#' @param a_lw Longwave absorptance
#' @param cp Heat capacity of dry air
#'
#' @return gR, leaf conductance to radiative heat transfer (m/sec)
#' @export
#'
#' @examples
conductance_radiative_heat <- function(Ta, Pa, rho, rho_mol, 
                                       a_lw=0.98, cp=1010, sb=5.67e-8) {
  2 * rho_mol * 4 * a_lw * sb * Ta^3 / (rho * cp)
}

#' Leaf boundary layer conductance to sensible heat transfer in
#' needle leaves.
#'
#' @param u Wind speed (m/sec)
#' @param d Characteristic leaf dimension (m)
#' @param rho_mol Molar air density (mol / m^3)
#'
#' @return gbH, conductance to sensible heat transfer (m/sec)
#' @export
#'
#' @examples
conductance_boundary_heat_needleleaf <- function(u, d, rho_mol) {
  0.006 * rho_mol * (u^0.6 / d^0.4)
}

#' Leaf boundary layer conductance to sensible heat transfer in
#' broad leaves.
#'
#' @param u Wind speed (m/sec)
#' @param d Characteristic leaf dimension (m)
#' @param rho_mol Molar air density (mol / m^3)
#'
#' @return gbH, conductance to sensible heat transfer (m/sec)
#' @export
#'
#' @examples
conductance_boundary_heat_broadleaf <- function(u, d, rho_mol) {
  0.0105 * rho_mol * (u^0.5 / d^0.5)
}

#' Psychrometric constant - relates partial pressure of water to
#' air temperature.
#'
#' @param Pa Air pressure (kPa)
#' @param cp Heat capacity of dry air (J/kg K) 
#' @param Mwr Ratio of molar mass of water to molar mass of dry air
#' @param l Heat of vaporization of water (J / mol)
#'
#' @return gamma, Psycrometric constant (kPa / K)
#' @export
#'
#' @examples
psychrometric_constant <- function(Pa, cp=1010, Mwr=0.622, l=40.8e3) {
  # Convert heat of vaporization to mass basis
  # 44.2e3 J/mol * 1/0.018 mol/kg = 2.45e6 J/kg
  l_mass <- l / 0.018
  cp * Pa / (l_mass * Mwr)
}

#' Slope of saturation vapor pressure with respect to temperature. Found by
#' differentiating Tetens's formula.
#'
#' @param Ta Air tempertaure (K)
#'
#' @return desat, Slope of esat w.r.t. air temperature (kPa / K)
#' @export
#'
#' @examples
esat_slope <- function(Ta) {
  a <- 17.27
  b <- 237.3
  Ta_c <- Ta - 273.15
  
  0.61078 * exp( (a * Ta_c) / (b + Ta_c) ) * a * b / (Ta_c+b)^2
}

#' Air-leaf temperature decoupling coefficient (Omega). Range is 0-1. When
#' omega is near 0, boundary layer conductance is high and leaf temperature
#' is well-coupled to air temperature. When omega is 1, boundary layer
#' conductance is low and latent heat flux is radiation-limited.
#'
#' @param gR Conductance to radiative heat transfer (mol / m2 sec)
#' @param gbH Boundary layer conductance to sensible heat (mol / m2 sec) 
#' @param gs Stomatal conductance (mol / m2 sec) 
#' @param d_esat Saturation VP vs. air temperature slope (kPa / K) 
#' @param gamma Psychrometric constant (kPa / K)
#' @param a_lw Longwave absorptance 
#'
#' @return omega, decoupling coefficient.
#' @export
#'
#' @examples
decoupling_coefficient <- function(gR, gbH, gs, d_esat, gamma, a_lw = 0.98) {
  eps <- d_esat / gamma
  
  omega <- (eps + 2 + (gR / gbH)) / (eps + ((gR+gbH) / gs) + gR/gbH)
  return(omega)
}

#' Equilibrium latent heat flux. This is the flux when leaf temperature is
#' decoupled from air temperature.
#'
#' @param Rn Net isothermal radiation (W / m2)
#' @param d_esat Slope of saturation VP vs. temperature (kPa / K)
#' @param gamma Psychrometric constant (kPa / K)
#' @param gR Leaf conductance to radiative heat transfer (mol / m2 sec)
#' @param gbH Leaf boundary layer conductance to heat transfer (mol / m2 sec) 
#'
#' @return LE_eq, equilibrium latent heat flux (W / m2)
#' @export
#'
#' @examples
latent_heat_equilibrium <- function(Rn, d_esat, gamma, gR, gbH) {
  eps <- d_esat / gamma
  le_eq <- (eps * Rn) / (eps + 1 + gR / gbH)
}

#' Total conductance. Assumes that boundary layer conductance to water vapor
#' is 1.08 times boundary layer heat conductance.
#'
#' @param gs Stomatal conductance (mol / m2 sec)
#' @param gbH Boundary layer heat conductance (mol / m2 sec) 
#'
#' @return gtot, total conductance (mol / m2 sec)
#' @export
#'
#' @examples
total_conductance <- function(gs, gbH) {
  (1/gs + 1/(gbH * 1.08)) ^ -1
}

#' Saturation vapor pressure from Tetens's formula.
#'
#' @param Ta Air temperature (K)
#'
#' @return esat, saturation vapor pressure (kPa)
#' @export
#'
#' @examples
saturation_vapor_pressure <- function(Ta) {
  Ta_c <- Ta - 273.15
  esat <- 0.61078 * exp((17.27 * Ta_c) / (Ta_c + 237.3))
  return(esat)
}

#' Vapor pressure deficit, the difference between the saturation vapor pressure
#' and the actual vapor pressure.
#'
#' @param Ta Air temperature (K)
#' @param RH Relative humidity (%, 0-100)
#'
#' @return vpd, vapor pressure deficit (kPa)
#' @export
#'
#' @examples
vapor_pressure_deficit <- function(Ta, RH) {
  vpd <- saturation_vapor_pressure(Ta) * (1 - RH/100)
  return(vpd)
}

#' Imposed latent heat flux. This is the latent heat flux when air temperature
#' and leaf temperature are well-coupled.
#'
#' @param gtot Total conductance (mol / m2 sec)
#' @param vpd Vapor pressure deficit (kPa) 
#' @param Pa Air pressure (kPa)
#' @param l Heat of vaporization of water (J / mol)
#'
#' @return LE_imp, stomatally-imposed latent heat flux.
#' @export
#'
#' @examples
latent_heat_imposed <- function(gtot, vpd, Pa, l=44.2e3) {
  le_imp <- gtot * l * vpd / Pa
  return(le_imp)
}

latent_heat <- function(omega, le_eq, le_imp) {
  le <- omega * le_eq + (1 - omega) * le_imp
  return(le)
}

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

# Useful constants
eb_constants <- list(
  cp=1010,    # Heat capacity of dry air J / kg K
  Mair=0.029, # Molar mass of dry air kg / mol
  a_lw=0.98,  # Longwave absorptance
  a_sw=0.50,  # Shortwave absorptance
  d=0.01,     # Characteristic dimension of leaves m
  MWr=0.622,  # Ratio of molar mass of water to dry air
  l=44.2e3    # Heat of vaporization of water J/mol *at 20 deg C*
)