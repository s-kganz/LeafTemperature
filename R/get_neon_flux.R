# Whole-canopy CO2 and H2O fluxes are not currently available through the
# Ameriflux API for all the sites we consider. NEON is currently working on
# getting data through the OneFlux pipeline, so we cannot generate a gap-filled
# dataset. However, we can pull the non gap-filled observations from NEON.

library(neonUtilities)
library(tidyverse)
library(foreach)

get_neon_flux <- function(site_meta, token, outdir) {
  flux_dpid <- "DP4.00200.001"
  
  # Determine site availability ----
  # We only use data in the non-provisional state since this means we can assume 
  # it has been QC'd for us.
  product_info <- getProductInfo(
    flux_dpid, token
  )
  
  product_avail <- product_info$siteCodes %>%
    filter(siteCode %in% site_meta$site_neon) %>%
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
  vars <- c(
    "timeBgn", "timeEnd",
    "data.fluxCo2.nsae.flux",
    "qfqm.fluxCo2.nsae.qfFinl",
    "data.fluxH2o.nsae.flux",
    "qfqm.fluxH2o.nsae.qfFinl"
  )
  
  flux_final <- foreach(i = 1:nrow(product_avail), .combine=rbind) %do% {
    
    this_site <- product_avail$siteCode[i]
    this_months <- product_avail$getMonths[[i]]
    # Drop site-months before 2022 for reproducibility
    before_2022 <- which(parse_number(str_split_i(this_months, "-", 1)) < 2022)
    this_months <- this_months[before_2022]

    cat("Now downloading", length(this_months),
        "months for site", this_site, "\n")
    
    foreach(m=this_months, .combine=rbind) %do% {

      suppressWarnings(zipsByProduct(
        dpID=flux_dpid,
        site=this_site,
        startdate=m,
        enddate=m,
        check.size=FALSE,
        token=token,
        savepath=tempdir()
      ))
      
      # Stack the table, select vars we care about
      fpath <- file.path(tempdir(), "filesToStack00200")
      flux_stack <- stackEddy(fpath, var=vars)
      
      flux_month <- flux_stack[[this_site]] %>%
        select(all_of(vars)) %>%
        mutate(site=this_site)
      
      # Delete detritus
      unlink(fpath, recursive=TRUE)
      stopifnot(!dir.exists(fpath))
      
      return(flux_month)
    }
  }
  
  write_if_not_exist(flux_final, file.path(outdir, "cross_site_flux.csv"))
}

