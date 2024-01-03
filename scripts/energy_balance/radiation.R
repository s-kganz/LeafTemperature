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