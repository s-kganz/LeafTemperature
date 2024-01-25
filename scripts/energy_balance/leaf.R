# Functions describing leaf physiology and transpiration for energy balance
# models.

source("scripts/energy_balance/constants.R")

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
                                       a_lw=0.98,
                                       cp=eb_constants_$cp,
                                       sb=eb_constants_$sb) {
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
  
  (eps + 2 + (gR / gbH)) / (eps + ((gR+gbH) / gs) + gR/gbH)
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
  (eps * Rn) / (eps + 1 + gR / gbH)
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
latent_heat_imposed <- function(gtot, vpd, Pa, l=eb_constants_$l) {
  gtot * l * vpd / Pa
}

#' Actual latent heat flux. Calculated as mixture of stomtally-imposed and
#' equilibrium endpoints. See Still et al. (2022) for details.
#'
#' @param omega Decoupling coefficient (0-1)
#' @param le_eq Equilibrium endpoint (W/m2)
#' @param le_imp Stomatally-imposed endpoint (W/m2)
#'
#' @return le, total latent heat flux (W/m2)
#' @export
#'
#' @examples
latent_heat <- function(omega, le_eq, le_imp) {
  omega * le_eq + (1 - omega) * le_imp
}

#' Sensible heat flux between a leaf and its boundary layer. Flux is multiplied
#' by 2 since both sides of leaf participate in heat transfer.
#'
#' @param gbH Leaf boundary layer conductance to heat transfer
#' @param Tl Leaf temperature (K)
#' @param Ta Air temperature (K)
#'
#' @return H, sensible heat flux (W/m2)
#' @export
#'
#' @examples
sensible_heat <- function(gbH, Tl, Ta) {
  2 * gbH * eb_constants_$cp * eb_constants_$Mair * (Tl - Ta)
}

#' Empirical model of stomatal conductance based on Medlyn et al. (2011).
#'
#' @param A Net photosynthesis (umol CO2 / m^2 s)
#' @param VPD Vapor pressure deficit (kPa)
#' @param CO2 Atmospheric CO2 concentration (umol CO2 / mol air == ppm) 
#' @param g1 Slope parameter (kPa^0.5, usually in range 1-10)
#' @param g0 Intercept parameter (mol / m^2 s, usually set to 0)
#'
#' @return
#' @export
#'
#' @examples
medlyn_gs <- function(A, VPD, CO2, g1, g0=0) {
  g0 + (1 + g1 / sqrt(VPD)) * (A / CO2)
}