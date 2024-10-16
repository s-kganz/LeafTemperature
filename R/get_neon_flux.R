# Whole-canopy CO2 and H2O fluxes are not currently available through the
# Ameriflux API for all the sites we consider. NEON is currently working on
# getting data through the OneFlux pipeline, so we cannot generate a gap-filled
# dataset. However, we can pull the non gap-filled observations from NEON.

#' Download eddy covariance data from NEON
#'
#' @param site_meta Site metadata table \code{\link{site_meta}}
#' @param token NEON API token. Can be NULL or NA if you don't have one.
#' @param outdir Directory to write stacked flux tables.
#' @rdname get_neon_flux
#' 
#' @description
#' 
#' There are two version of this function: `get_neon_flux_long()` and `get_neon_flux_simple()`. Due to cross-platform issues with file manipulation, the final flux table lives in a Google Cloud bucket. `get_neon_flux_long()` will try to download everything to your machine for the sake of reproducibility. `get_neon_flux_simple()` just downloads the final CSV. ¯\_(ツ)_/¯
#' 
#' This function downloads eddy covariance data from NEON, stacks the level 4 NSAE data and QC flags, and then writes that to a CSV in `outdir`. For growing season data only, this downloads about 60 GB total, but uses <500 MB on disk at any one time by deleting intermediate files. This is suboptimal, but is currently recommended for large flux downloads (see \href{https://github.com/NEONScience/NEON-utilities/issues/131}{this github issue}).
#' 
#' The function will try to delete intermediate folders generated during the stacking process. But, deleting files does not work cons
#' 
#' For reproducibility with the manuscript, this script also excludes any flux products after 2021.
#' 
#'
#' @export
#'
get_neon_flux_long <- function(site_meta, token, outdir) {
  flux_dpid <- "DP4.00200.001"
  
  # Determine site availability ----
  # We only use data in the non-provisional state since this means we can assume 
  # it has been QC'd for us.
  product_info <- neonUtilities::getProductInfo(
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
  
  temp_flux_loc <- "temp_flux_loc"
  if (temp_flux_loc != tempdir()) {
    dir.create(temp_flux_loc)
  }
  
  flux_final <- foreach(i = 1:nrow(product_avail), .combine=rbind) %do% {
    
    this_site <- product_avail$siteCode[i]
    this_months <- product_avail$getMonths[[i]]
    # Drop site-months before 2022 for reproducibility
    before_2022 <- which(parse_number(str_split_i(this_months, "-", 1)) <= 2022)
    this_months <- this_months[before_2022]

    message(paste(
      "Now downloading",
      length(this_months),
      "months for site",
      this_site
    ))
    
    foreach(m=this_months, .combine=rbind) %do% {

      suppressWarnings(neonUtilities::zipsByProduct(
        dpID=flux_dpid,
        site=this_site,
        startdate=m,
        enddate=m,
        check.size=FALSE,
        token=token,
        savepath=temp_flux_loc
      ))
      
      # Stack the table, select vars we care about
      fpath <- file.path(temp_flux_loc, "filesToStack00200")
      flux_stack <- neonUtilities::stackEddy(fpath, var=vars)
      
      flux_month <- flux_stack[[this_site]] %>%
        select(all_of(vars)) %>%
        mutate(site=this_site)
      
      # Delete detritus
      rm(flux_stack)
      unlink(fpath, recursive=TRUE, force=TRUE)
      if (dir.exists(fpath)) {
        message("Failed to delete")
        stop()
      }

      return(flux_month)
    }
  }
  
  write_if_not_exist(flux_final, file.path(outdir, "cross_site_flux.csv"))
}

#' @rdname get_neon_flux
#' @export
get_neon_flux_simple <- function(outdir) {
  url <- "https://storage.googleapis.com/neon-flux-table/cross_site_flux.csv"
  fname <- basename(url)
  
  download.file(url, file.path(outdir, fname))
}
