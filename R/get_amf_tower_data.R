# Download tower micrometeorology data from Ameriflux. Also, collect all of
# the radiometer canopy temperature data in one table.

#' Download tower data from Ameriflux
#'
#' @param site_meta Site metadata table \code{\link{site_meta}}.
#' @param outdir Directory to write tower CSVs. A subfolder named `amf_tower_data` will be created here.
#' @param user Ameriflux username.
#' @param email Ameriflux email.
#'
#' @description
#' This function downloads the current Ameriflux BASE product for each site listed in `site_meta`. It also isolates all `T_CANOPY` observations and writes these to a separate table. For the 11 sites listed in \code{\link{site_meta}}, about 2GB of storage is needed.
#' 
#'
#' @export
get_amf_tower_data <- function(site_meta, outdir, user, email) {
  if (!dir.exists(outdir)) dir.create(outdir)
  dir.create(file.path(outdir, "amf_tower_data"))
  
  flux_files <- amerifluxr::amf_download_base(
    user,
    email,
    site=site_meta$site_ameriflux,
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
    this_i    <- which(this_site == site_meta$site_ameriflux)
    this_utc  <- site_meta$utc_offset[this_i]
    
    message(paste(f, "\n", this_site, this_i, this_utc))
    
    this_f <- amerifluxr::amf_read_base(
      file = f,
      unzip = TRUE,
      parse_timestamp = TRUE
    ) %>% 
      amerifluxr::amf_filter_base() %>%
      # The timestamp is recorded in local standard time *but is parsed as if it is UTC*
      # so we have to manually convert to UTC
      # https://github.com/chuhousen/amerifluxr/issues/95
      mutate(
        TIMESTAMP = TIMESTAMP - minutes(15) - hours(this_utc)
      ) %>%
      # Drop extra timekeeping cols
      select(-TIMESTAMP_START, -TIMESTAMP_END, -YEAR, -MONTH, -DAY, -DOY, -HOUR, -MINUTE) %>%
      # Filter observations to prior to 2022 for reproducibility
      filter(year(TIMESTAMP) <= 2022) %>%
      write_if_not_exist(
        file.path(outdir, "amf_tower_data", paste0(this_site, ".csv"))
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
    write_if_not_exist(file.path(outdir, "cross_site_tcan.csv"))
  
}

