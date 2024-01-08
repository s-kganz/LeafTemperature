# Whole-canopy CO2 and H2O fluxes are not currently available through the
# Ameriflux API for all the sites we consider. NEON is currently working on
# getting data through the OneFlux pipeline, so we cannot generate a gap-filled
# dataset. However, we can pull the non gap-filled observations from NEON.

library(neonUtilities)
library(tidyverse)
library(foreach)

# Determine site availability ----
# We only use data in the non-provisional state since this means we can assume 
# it has been QC'd for us.

site_info <- read_csv("data_working/neon_site_metadata.csv")
token <- readChar("data_in/neon_token", file.info("data_in/neon_token")$size)
flux_dpid <- "DP4.00200.001"

product_info <- getProductInfo(
  flux_dpid, token
)

product_avail <- product_info$siteCodes %>%
  filter(siteCode %in% site_info$site_neon) %>%
  select(-availableMonths) %>%
  unnest(availableReleases) %>%
  filter(release != "PROVISIONAL") %>%
  # These data can be quite large so only grab the growing season months
  mutate(
    getMonths = map(
      availableMonths,
      function(v) {
        v[which(str_detect(v, "-0[5-9]"))]
      }
    )
  )

# Download monthly flux data ----

out <- tempdir()

n_months <- 0
for (i in 1:nrow(product_avail)) {
  this_site <- product_avail$siteCode[i]
  this_months <- product_avail$getMonths[[i]]
  n_months <- n_months + length(this_months)
  
  cat("Now downloading", length(this_months), 
      "months for site", this_site, "\n")
  
  vapply(this_months, function(m) {
    cat("\t", m, "\n")
    suppressWarnings(zipsByProduct(
      dpID=flux_dpid,
      site=this_site,
      startdate=m,
      enddate=m,
      check.size=FALSE,
      token=token,
      savepath=out
    ))
    # Throw away return value to make vapply hapy
    return(0)
  }, 0)
}

# Stack tables ----

fpath <- file.path(out, "filesToStack00200", fsep="\\")

n_files <- length(list.files(fpath, pattern="*.zip"))

stopifnot(n_files == n_months)

vars <- c(
  "timeBgn", "timeEnd",
  "data.fluxCo2.nsae.flux",
  "qfqm.fluxCo2.nsae.qfFinl",
  "data.fluxH2o.nsae.flux",
  "qfqm.fluxH2o.nsae.qfFinl"
)

flux_stack <- stackEddy(fpath, var=vars)

flux_clean <- flux_stack[product_avail$siteCode]

cat("# Observations per site:\n")

flux_final <- foreach(site=names(flux_clean), data=flux_clean, 
                      .combine=rbind) %do% 
  {
    data_clean <- data %>%
      select(all_of(vars)) %>%
      mutate(site=site) %>%
      filter(qfqm.fluxH2o.nsae.qfFinl == 0,
             qfqm.fluxCo2.nsae.qfFinl == 0,
             !is.nan(data.fluxCo2.nsae.flux),
             !is.nan(data.fluxH2o.nsae.flux))
    
    cat("\t", site, "\t", nrow(data_clean), "\n")
    
    return(data_clean)
  }

write_csv(flux_final, "data_out/cross_site_flux_qc.csv")
