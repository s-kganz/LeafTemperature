
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