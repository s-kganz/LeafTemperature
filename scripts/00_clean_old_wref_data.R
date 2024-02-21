library(tidyverse)
source("scripts/tower_util.R")
source("scripts/energy_balance/air.R")

wref <- readxl::read_xlsx(
  "data_in/chris_still_wref_tir/WR_TIRtemp,Met,Flux_1hr_140101-151231_160425-with-metadata-column-headers.xlsx",
  sheet=2,
  skip=1
) %>%
  # Time zone accounting
  mutate(TIMESTAMP = parse_date_time(paste(Yr, DOY, Hr), "yjh", tz="Etc/GMT-8"),
         TIMESTAMP = with_tz(TIMESTAMP, "UTC"))

# Clean up names and units to match Ameriflux as closely as possible
wref_clean <- wref %>%
  select(TIMESTAMP, matches("Tcan_"), everything(), -matches("Ts_\\d"),
         -matches("S1_\\d"), -Yr, -Mon, -DOY, -Hr, -tDOY, -`dir_02(o)`,
         -matches("_2nd"), -TOTAL, -n_Tot, -factor_CO2, -factor_H2O,
         -`H2O_70(g m-3)`, -matches("par_\\d"), -`F_H2O_70(g m-2 s-1)`,
         -`dir_70(o)`) %>%
  rename(
    P = "prcpt1 (mm)",
    PA = "p_70(kPa)",
    TA_02 = "Ta_02(oC)",
    RH_02 = "RH_02(%)",
    WS_02 = "u_02(m s-1)",
    G = "G(W m-2)",
    SW_IN = "S_in(W m-2)",
    SW_OUT = "S_out(W m-2)",
    LW_IN = "L_in(W m-2)",
    LW_OUT = "L_out(W m-2)",
    NEE = "F_CO2_70(mg m-2 s-1)",
    H = "F_Hs_70(W m-2)",
    LE = "F_LE_70(W m-2)",
    USTAR = "u_star(m s-1)",
    TA_70 = "Ta_70(oC)",
    RH_70 = "RH_70(%)",
    WS_70 = "u_70(m s-1)",
    SW_DIF = "DIFF",
    CO2 = "CO2_70(mg m-3)",
    Tcan_avg = "Tcan_AvgS0-9"
  ) %>%
  mutate(
    # mg CO2 --> umol CO2
    NEE = NEE * 1e3 / 44.01,
    CO2_umol_m3 = CO2 * 1e3 / 44.01,
    rho_mol = dry_air_molar_density(TA_70+273.15, PA),
    CO2_ppm = CO2_umol_m3 / rho_mol
    #rho_mol = dry_air_molar_density(TA_70+273.15, )
  ) %>%
  select(-CO2_umol_m3, -rho_mol, -CO2) %>%
  rename(CO2 = CO2_ppm)

write_if_not_exist(wref_clean, "data_out/old_wref_clean.csv")
