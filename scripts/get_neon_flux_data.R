# Download all the flux data needed for the multisite synthesis
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
  this_f <- amf_read_base(
    file = f,
    unzip = TRUE,
    parse_timestamp = TRUE
  ) %>% 
    amf_filter_base() %>%
    write_csv(
      file.path(outdir, paste0(str_split_i(basename(f), "_", 2), ".csv"))
    )
}
