# Within-canopy absorbed radiation models from Wang and Leuning (1998)

source("scripts/energy_balance/constants.R")

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
net_isothermal_radiation <- function(SW_abs, LW_abs_dn, Ta, a_lw) {
  sb <- eb_constants_$sb
  # Multiply by 2 because both sides of the leaf radiate heat.
  Rn <- SW_abs + LW_abs_dn - outgoing_longwave_radiation(Ta, ef=a_lw)
}

#' Net radiation, but allowing for nonisothermal conditions.
#'
#' @param SW_abs Absorbed downwelling shortwave radiation (W/m2)
#' @param LW_abs_dn Absorbed downwelling longwave radiation (W/m2)
#' @param Tl Leaf temperature (K)
#' @param a_sw Leaf absorptance to shortwave radiation (0-1)
#' @param a_lw Leaf absorptance to longwave radiation (0-1)
#'
#' @return Rnet, net radiation (W/m2)
#' @export
#'
#' @examples
net_nonisothermal_radiation <- function(SW_abs, LW_abs_dn, Tl, a_lw=0.98) {
  SW_abs + LW_abs_dn - outgoing_longwave_radiation(Tl, ef=a_lw)
}

#' Absorbed shortwave radiation within a tree canopy. Comes from Wang and
#' Leuning (1998) Eqs. B2-3
#'
#' @param SW_beam Top-of-canopy direct shortwave radiation (W/m2)
#' @param SW_dif Top-of-canopy diffuse shortwave radiation (W/m2)
#' @param l Cumulative leaf area index from top of canopy (m2 / m2)
#' @param kb Attenuation coefficient for black leaves under direct radiation
#' @param kd Attenuation coefficient for black leaves under diffuse radiation
#' @param kb_real Empirical attenuation coefficient under direct radiation
#' @param kd_real Empirical attenuation coefficient under diffuse radiation
#' @param s Leaf scattering coefficient = 1 - absorptance (0-1)
#' @param rho_b Canopy reflectance coefficient for direct radiation
#' @param rho_d Canopy reflectance coefficient for diffuse radiation
#' @param aslist Return partitioned sun/shade fluxes as a list?
#'
#' @return q_total, total absorbed shortwave radiation of sunlit and 
#' shaded leaves (W/m2)
#' @export
#'
#' @examples
absorbed_shortwave_radiation <- function(SW_beam, SW_dif, l,
                                         kb_real, kd_real, kb, kd=0.8,
                                         s=0.5, rho_b=0.1, rho_d=0.1, 
                                         aslist=FALSE) {
  q_shade <- (
    SW_dif * kd_real * (1 - rho_d) * exp(-kd_real * l) +
    SW_beam * (kb_real * (1-rho_b) * exp(-kb_real * l)) -
    kb * (1 - s) * exp(-kb * l)
  )
  
  q_sun <- q_shade + SW_beam * kb * (1 - s)
  
  f_sun <- exp(-kb * l)
  
  q_total <- q_sun * f_sun + q_shade * (1 - f_sun)
  
  if (aslist) {
    return(list(
      q_sun=q_sun,
      q_shade=q_shade,
      f_sun=f_sun,
      q_total=q_total
    ))
  } else {
    return(q_total)
  }
}

#' Absorbed shortwave radiation within a tree canopy. Comes from Wang and
#' Leuning (1998) Eqs. B8-11
#'
#' @param Ta Air temperature (K)
#' @param l Cumulative leaf area index from top of canopy (m2 / m2)
#' @param L Total leaf area index for site
#' @param ef Emissivity of foliage
#' @param ea Effective emissivity of the sky
#' @param es Emissivity of the soil
#' @param kb Attenuation coefficient for black leaves under beam radiation
#' @param kd Attenuation coefficient for black leaves under diffuse radiation
#' @param aslist Return partitioned sun/shade fluxes as a list?
#'
#' @return q_total, total absorbed longwave radiation of sunlit and 
#' shaded leaves (W/m2)
#' @export
#'
#' @examples
absorbed_longwave_radiation <- function(Ta, l, L,
                                        kb, ea, ef=0.95, es=0.9,
                                        kd=0.8, aslist=FALSE) {
  # Effective emissivity of downwelling radiation
  e_down <- (
    ef -
    (ef - ea) * exp(-kd * l) +
    (1 - ef) * (1 - exp(-kd * l)) * (ef - (ef - es) * exp(-kd * (L - l)))
  )
  
  # Effective emissivity of upwelling radiation
  e_up <- (
    ef - 
    (ef - es) * exp(-kd * (L - l)) +
    (ef - (ef - ea) * exp(-kd * l)) * ((1 - ef) * (1 - exp(-kd * (L - l)))) +
    (1 - es) * exp(-2 * kd * (L - l)) +
    ef * (1 - es) * exp(-kd * (L - l)) * (1 - exp(-kd * (L - l)))
  )
  
  # Absorbed radiation
  sigma <- eb_constants_$sb # SB const
  q_up <- e_up * sigma * (Ta^4)
  q_down <- e_down * sigma * (Ta^4)
  q_net <- q_up - q_down
  
  # Partition into sun and shade
  f_sun <- exp(-kb * l)
  
  q_sun <- f_sun * q_net
  q_shade <- (1-f_sun) * q_net
  
  if (aslist) {
    return(list(
      e_down=e_down,
      e_up=e_up,
      q_sun=q_sun,
      q_shade=q_shade,
      f_sun=f_sun,
      q_net=q_net
    ))
  } else {
    return(q_net)
  }
}

#' Downwelling longwave radiation. This essentially takes the same approach
#' as the version from Wang and Leuning implemented above but only calculates
#' the downwelling component.
#'
#' @param Ta Air temperature (K)
#' @param l Cumulative leaf area index from top of canopy (m2 / m2)
#' @param L Total leaf area index for site
#' @param ef Emissivity of foliage
#' @param ea Effective emissivity of the sky
#' @param es Emissivity of the soil
#' @param kb Attenuation coefficient for black leaves under beam radiation
#' @param kd Attenuation coefficient for black leaves under diffuse radiation
#'
#' @return L_dn, downwelling longwave radiation
#' @export
#'
#' @examples
absorbed_downwelling_longwave_radiation <- function(Ta, l, L,
                                           kb, ea, ef=0.95, es=0.9,
                                           kd=0.8, aslist=FALSE) {
  # Effective emissivity of downwelling radiation
  e_down <- (
    ef -
      (ef - ea) * exp(-kd * l) +
      (1 - ef) * (1 - exp(-kd * l)) * (ef - (ef - es) * exp(-kd * (L - l)))
  )
  
  # Multiplying by 2 here because longwave radiation is diffuse and radiated
  # in all directions. So, the underside of leaves also receives this flux.
  LW_abs_dn <- 2 * eb_constants_$sb * e_down * Ta^4
  
  if (aslist) {
    return(list(
      e_down=e_down,
      LW_abs_dn=LW_abs_dn
    ))
  } else {
    return(LW_abs_dn)
  }
}

#' Outgoing (upwelling) longwave radiation from a leaf. This version just models
#' the leaf as a blackbody since the upwelling component from Wang and Leuning
#' (1998) is essentially equivalent to this. The flux is doubled because both
#' sides of the leaf participate in heat transfer.
#'
#' @param Tl Leaf temperature (K)
#' @param ef Effective foliage emissivity (0-1)
#'
#' @return L_up, upwelling longwave radiation
#' @export
#'
#' @examples
outgoing_longwave_radiation <- function(Tl, ef=0.95) {
  2 * ef * eb_constants_$sb * Tl^4
}

if (sys.nframe() == 0) {
  # Sensible defaults
  L <- 5.5
  l <- seq(0, L, by=0.1)
  
  ef <- 0.95
  es <- 0.9
  ea <- 0.83
  kd <- 0.8
  kb <- 0.5
  kd_real <- 0.8
  kb_real <- 0.7
  
  SW_beam <- 1
  SW_dif <- 0.5
  
  q_sw <- absorbed_shortwave_radiation(SW_beam, SW_dif, l, kb_real, kd_real,
                                       kb, kd, s=0.5, rho_b=0.04, rho_d=0.06)
  
  q_lw <- absorbed_longwave_radiation(298, l, L, kb, ea)
  
  plot(l, q_sw)
  plot(l, q_lw)
}
