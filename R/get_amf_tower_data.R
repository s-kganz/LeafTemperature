# Download tower micrometeorology data from Ameriflux. Also, collect all of
# the radiometer canopy temperature data in one table.
library(tidyverse)
library(amerifluxr)

get_amf_tower_data <- function(sites, outdir, user, email) {
  if (!dir.exists(outdir)) dir.create(outdir)
  
  flux_files <- amf_download_base(
    user,
    email,
    site=sites$site_ameriflux,
    data_policy="CCBY4.0",
    agree_policy=TRUE,
    intended_use = "other",
    intended_use_text = "canopy temperature synthesis",
    verbose = TRUE,
    out_dir = tempdir()
  )
  
  amf_tcan <- data.frame(TIMESTAMP=c(), SITE=c(), TCAN_VAR=c(), TCAN=c())
  
  for (f in flux_files) {
    this_site <- str_split_i(basename(f), "_", 2)
    this_i    <- which(this_site == sites$site_ameriflux)
    this_utc  <- sites$utc_offset[this_i]
    
    cat(f, "\n", this_site, this_i, this_utc, "\n")
    
    this_f <- amf_read_base(
      file = f,
      unzip = TRUE,
      parse_timestamp = TRUE
    ) %>% 
      amf_filter_base() %>%
      # The timestamp is recorded in local standard time *but is parsed as if it is UTC*
      # so we have to manually convert to UTC
      # https://github.com/chuhousen/amerifluxr/issues/95
      mutate(
        TIMESTAMP = TIMESTAMP - minutes(15) - hours(this_utc)
      ) %>%
      # Drop extra timekeeping cols
      select(-TIMESTAMP_START, -TIMESTAMP_END, -YEAR, -MONTH, -DAY, -DOY, -HOUR, -MINUTE) %>%
      write_if_not_exist(
        file.path(outdir, paste0(this_site, ".csv"))
      )
    
    # Collect all the Tcan values and add to the table
    this_tcan <- this_f %>%
      select(TIMESTAMP, matches("T_CANOPY_\\d_\\d_\\d")) %>%
      mutate(SITE=this_site) %>%
      pivot_longer(matches("T_CANOPY_\\d_\\d_\\d"), names_to="TCAN_VAR",
                   values_to="TCAN")
    
    amf_tcan <- bind_rows(amf_tcan, this_tcan)
  }
  
  amf_tcan %>%
    drop_na() %>%
    filter(month(TIMESTAMP) %in% 5:9) %>%
    write_if_not_exist("data_out/cross_site_tcan.csv")
  
}

