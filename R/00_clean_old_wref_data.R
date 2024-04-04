# This script was written for use with the flux and TIR data from the Still 
# et al. PNAS paper. The source data for this is available at
# https://zenodo.org/records/6862565. The script is here in case anyone else
# wants to play with these data, but the analysis is not included in our paper.

library(tidyverse)
library(readxl)
source("scripts/tower_util.R")
source("scripts/energy_balance/air.R")

wref <- read_csv(
  "data_in/chris_still_wref_tir/still_WR_flux_data.csv",
  col_types=cols(CO2_1=col_double(), CO2_2=col_double())
) %>%
  # Time zone accounting. For some reason we have to flip the sign of the UTC
  # offset. This is apparently a thing in core R, idk why it is like that.
  # https://stackoverflow.com/questions/57518552/problem-with-converting-unix-time-zone-in-lubridate
  mutate(
    TIMESTAMP_START = as.character(TIMESTAMP_START),
    TIMESTAMP_START = ymd_hm(TIMESTAMP_START, tz="Etc/GMT+8"),
  )

# Verify timezone is correct.
# wref %>%
#   mutate(tz_local = with_tz(TIMESTAMP_START, tz="US/Pacific")) %>%
#   ggplot(aes(x=hour(tz_local), y=SW_IN)) + geom_boxplot(aes(group=hour(tz_local)))

# Clean up names and units to match Ameriflux as closely as possible
wref_clean <- wref %>%
  rename(TIMESTAMP = TIMESTAMP_START, H = H.x, LE = LE.x) %>%
  select(TIMESTAMP, matches("Tcan_S\\d"), everything(), -FC, -SC, -TS_1, -TS_2,
         -SLE, -SH, -WD, -SWC_1, -SWC_2, -FAPAR, -ZL, -DateTime, -DOY.x, -(Y:s),
         -LE.y, -H.y, -Tair, -Tsoil, -Reco_DT, -GPP_DT, -Yr, -Mon, -DOY.y,
         -Hr, -tDOY, -TIMESTAMP_END, -Tcan_Avg, -SkyT, -Tcan_Avg_corr, -date,
         -DOY, -hour, -rH, -VPD, -APAR, -PPFD_DIF, -H2O, -RECO_PI, -VPD_PI,
         -NEE_PI, -Rg)

# CO2 is totally NA so we have to pull that from the Excel sheet
wref_xl <- readxl::read_xlsx(
  "data_in/chris_still_wref_tir/WR_TIRtemp,Met,Flux_1hr_140101-151231_160425-with-metadata-column-headers.xlsx",
  sheet=2,
  skip=1,
  col_types="numeric"
) %>%
  # Time zone accounting. For some reason we have to flip the sign of the UTC
  # offset. This is apparently a thing in core R, idk why it is like that.
  # https://stackoverflow.com/questions/57518552/problem-with-converting-unix-time-zone-in-lubridate
  mutate(TIMESTAMP = parse_date_time(paste(Yr, DOY, Hr), "yjh", tz="Etc/GMT+8"),
         TIMESTAMP = with_tz(TIMESTAMP, "UTC")) %>%
  select(TIMESTAMP, `CO2_70(mg m-3)`) %>%
  rename(CO2_mg_m3 = `CO2_70(mg m-3)`) %>%
  drop_na()

wref_clean_co2 <- wref_clean %>%
  inner_join(wref_xl, by="TIMESTAMP") %>%
  mutate(
    # mg CO2 --> umol CO2
    CO2_umol_m3 = CO2_mg_m3 * 1e3 / 44.01,
    rho_mol = dry_air_molar_density(TA+273.15, PA),
    CO2_ppm = CO2_umol_m3 / rho_mol
    #rho_mol = dry_air_molar_density(TA_70+273.15, )
  ) %>%
  select(-CO2_umol_m3, -rho_mol, -CO2_mg_m3, -CO2_1, -CO2_2) %>%
  rename(CO2 = CO2_ppm)

write_if_not_exist(wref_clean_co2, "data_out/old_wref_clean.csv")
