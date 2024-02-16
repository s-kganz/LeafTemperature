library(tidyverse)
library(REddyProc)
# For VPD function
source("scripts/energy_balance/air.R")
source("scripts/tower_util.R")

process_site_year <- function(site, lat, lon, flux, tower) {
  flux_tower <- flux %>%
    select(timeBgn, data.fluxCo2.nsae.flux) %>%
    left_join(
      tower %>%
        select(TIMESTAMP, SW_IN_1_1_1, TA_1_1_1, RH, USTAR),
      by=c("timeBgn"="TIMESTAMP")
    ) %>%
    setNames(c("DateTime", "NEE", "Rg", "Tair", "RH", "Ustar")) %>%
    # REddyProc expects VPD in hPa instead of kPa like everyone else >:(
    mutate(VPD = vapor_pressure_deficit(Tair + 273.15, RH) * 10,
           # We don't apply the NEON QC flags until after gap-filling; these
           # bounds ensure we don't use really bad data.
           Rg = ifelse(Rg < 0, 0, Rg),
           VPD = ifelse(VPD < 0, 0, VPD),
           NEE = ifelse(NEE < -50, NA, NEE),
           NEE = ifelse(NEE > 100, NA, NEE))
  
  flux_tower <- filterLongRuns(flux_tower, "NEE")
  
  EProc <- sEddyProc$new(
    site, flux_tower, c('NEE','Rg','Tair','VPD', 'Ustar')
  )
  
  EProc$sEstimateUstarScenarios(
    nSample = 100L, probs = c(0.05, 0.5, 0.95)
  )
  
  # Gap fill where Ustar is below the cutoff
  EProc$sMDSGapFillUStarScens("NEE")
  
  # Reco and GPP partitioning
  # Timestamps are assumed to be in UTC
  EProc$sSetLocationInfo(LatDeg = lat, LongDeg = lon, TimeZoneHour = 0)  
  EProc$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)     
  EProc$sMDSGapFill('VPD', FillAll = FALSE,  minNWarnRunLength = NA)     
  EProc$sFillVPDFromDew() # fill longer gaps still present in VPD_f
  
  EProc$sMRFluxPartitionUStarScens(parsE0Regression=list(TempRange=3))
  

  result <- EProc$sExportResults()
  
  if ("R_ref_U50" %in% names(result)) {
    # Partitioning succeeded!
    result %>%
      # Collect variables of interest
      select(Ustar_U50_Thres, R_ref_U50, E_0_U50, Reco_U50, GPP_U50_f) %>%
      # Join these back with the original table. This assumes REddyProc doesn't
      # eat any rows.
      cbind(flux, .) %>%
      return()
  } else {
    # Partitioning failed :(
    # Just fill with NAs so the script will still finish.
    result %>%
      select(Ustar_U50_Thres) %>%
      mutate(R_ref_U50=NA,
             E_0_U50=NA,
             Reco_U50=NA,
             GPP_U50_f=NA) %>%
      cbind(flux, .) %>%
      return()
  }
}

get_tower_data <- function(site, year) {
  data_dir <- "data_working/neon_flux"
  file <- file.path(data_dir, str_c(site, ".csv"))
  
  read_csv(file, col_types=cols()) %>%
    filter(year(TIMESTAMP) == year) %>%
    select(TIMESTAMP, SW_IN_1_1_1, TA_1_1_1, RH, USTAR) %>%
    mutate(site=site) %>%
    return()
}

if (sys.nframe() == 0) {
  site_meta <- read_csv("data_working/neon_site_metadata.csv")
  all_flux <- read_csv("data_out/cross_site_flux.csv") %>%
    # N.b. can't drop any rows here otherwise REddyProc complains about non-
    # equidistant time steps
    select(-timeEnd) %>%
    mutate(year = year(timeBgn)) %>%
    group_by(site, year) %>%
    group_split()
  
  processed_flux <- lapply(all_flux, function(flux) {
    # Check that there are at least 90 days of data. If not, REddyProc won't
    # process it, so we have to throw it out
    
    # Get site/year
    site_neon <- flux$site[1]
    year <- flux$year[1]
    print(c(site_neon, year))
    
    timediff <- max(flux$timeBgn) - min(flux$timeBgn)
    if (timediff < days(90)) {
      warning(site_neon, "(", year, ")", "has less than 90 days of data. ",
              "Returning NA.")
      return(NA)
    }
    
    if (sum(!is.na(flux$data.fluxCo2.nsae.flux)) < 3000) {
      warning(site_neon, " (", year, ") ", "has insufficient non-NA data. ",
                  "Returning NA.")
      return(NA)
    }
    
    # Lookup other params from metadata table
    meta_idx <- match(site_neon, site_meta$site_neon)
    site_amf <- site_meta$site_ameriflux[meta_idx]
    site_lat <- site_meta$tower_lat[meta_idx]
    site_lon <- site_meta$tower_lon[meta_idx]
    
    # Pull tower data
    tower <- get_tower_data(site_amf, year)
    
    # Pass to REddyProc pipeline
    process_site_year(site_amf, site_lat, site_lon, flux, tower)
  })
  
  # Drop NA results and bind together in one big table
  processed_flux_df <- processed_flux[lapply(processed_flux, class) == "data.frame"] %>%
    bind_rows()
  
  processed_flux_df_qc <- processed_flux_df %>%
    # Select only high-quality flux measurements
    filter(qfqm.fluxCo2.nsae.qfFinl == 0, qfqm.fluxH2o.nsae.qfFinl == 0) %>%
    # Drop years where partitioning failed
    filter(!is.na(R_ref_U50)) %>%
    # Drop observations where modeled partitioning does not match measured NEE
    mutate(NEE_model = Reco_U50 - GPP_U50_f,
           reldiff = (NEE_model - data.fluxCo2.nsae.flux) / data.fluxCo2.nsae.flux) %>%
    filter(abs(reldiff) < 0.1) %>%
    # Drop observations where GPP < 0 since this is nonphysical
    filter(GPP_U50_f >= 0) %>%
    select(-NEE_model, -reldiff)
  
  processed_flux_df_qc %>%
    pull(site) %>%
    table()
  
  processed_flux_df_qc %>%
    select(-year) %>%
    write_if_not_exist("data_out/cross_site_flux_partition_qc.csv")
}
