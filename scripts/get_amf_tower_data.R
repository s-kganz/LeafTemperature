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

for (i in 1:length(flux_files)) {
  this_f <- amf_read_base(
    file = flux_files[i],
    unzip = TRUE,
    parse_timestamp = FALSE
  ) %>% 
    amf_filter_base() %>%
    # The timestamp starts in local standard time, convert to UTC so we can
    # save the file without anything changing.
    mutate(
      TIMESTAMP = ymd_hm(TIMESTAMP_START) - hours(sites$utc_offset[i]) + minutes(15)
    ) %>%
    select(-TIMESTAMP_START, TIMESTAMP_END) %>%
    write_csv(
      file.path(outdir, paste0(str_split_i(basename(flux_files[i]), "_", 2), ".csv"))
    )
}
