# Physical constants for energy balance models

if (!exists("eb_constants_")) {
  eb_constants_ <- list(
    cp=1010,    # Heat capacity of dry air J / kg K
    Mair=0.029, # Molar mass of dry air kg / mol
    MWr=0.622,  # Ratio of molar mass of water to dry air
    l=44.2e3,   # Heat of vaporization of water J/mol *at 20 deg C*
    sb=5.67e-8, # Stefan-Boltzmann constant W m-2 K-4
    R=8.314     # Universal gas constant J / mol K
  )
}