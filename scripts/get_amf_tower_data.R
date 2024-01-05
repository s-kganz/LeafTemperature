# Download tower micrometeorology data from Ameriflux
library(tidyverse)
library(amerifluxr)

sites <- read_csv("data_working/neon_site_metadata.csv")
outdir <- file.path("data_working", "neon_flux")

flux_files <- amf_download_base(
  "ganzk",
  "ganzk@uw.edu",
  site=sites$site_ameriflux,
  data_policy="CCBY4.0",
  agree_policy=TRUE,
  intended_use = "other",
  intended_use_text = "canopy temperature synthesis",
  verbose = TRUE,
  out_dir = tempdir()
)

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
    write_csv(
      file.path(outdir, paste0(this_site, ".csv"))
    )
}
