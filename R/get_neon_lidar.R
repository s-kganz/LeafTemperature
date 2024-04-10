#' Download NEON aerial LiDAR near flux towers
#'
#' @param site_meta Site metadata table \code{\link{site_meta}}.
#' @param outdir Output directory where processed point clouds are saved.
#'
#' @description
#' This function downloads all aerial LiDAR tiles within 50 m of the coordinates in `site_meta`, normalizes height, and then clips the points together into one 50 m circle. 
#' 
#'
#' @export
get_neon_lidar <- function(site_meta, outdir) {
  if (!dir.exists(outdir)) dir.create(outdir)

  product <- "DP1.30003.001"
  lidar_avail <- neonUtilities::getProductInfo(product)
  
  buffer_radius <- 50 # m

  for(i in 1:nrow(site_meta)) {
    site_id <- site_meta$site_neon[i]
    lat <- site_meta$tower_lat[i]
    lon <- site_meta$tower_lon[i]
    utm <- site_meta$utm_zone[i]
    
    # Find the most recent acquisition
    site_index <- match(site_id, lidar_avail$siteCodes$siteCode)
    avail_releases <- lidar_avail$siteCodes$availableReleases[[site_index]] %>%
      filter(release != "PROVISIONAL")
    acq_year_mo <- tail(avail_releases$availableMonths[[1]], n=1)
    acq_year <- strsplit(acq_year_mo, "-")[[1]][1]
    
    # Convert lat lon to UTM
    utm_proj4 <- paste("+proj=utm ", "+zone=", utm, " +datum=NAD83 ", "+units=m",
                       sep="")
    utm_coord <- terra::project(
      matrix(c(lon, lat), nrow=1), 
      from=terra::crs("epsg:4326"), to=terra::crs(utm_proj4)
    )
    
    utm_x <- utm_coord[1, 1]
    utm_y <- utm_coord[1, 2]
    
    # Download files to a tempdir
    tdir <- tempdir()
    
    neonUtilities::byTileAOP(
      product,
      site_id,
      acq_year,
      utm_x, # easting
      utm_y, # northing
      buffer=buffer_radius,
      check.size=FALSE,
      savepath=tdir
    )
    
    # Make a las catalog from the files
    # print("Reading downloaded LAS files...")
    lascat <- lidR::readLAScatalog(
      list.files(file.path(tdir, product), pattern="*.laz", recursive=TRUE,
                 full.names=TRUE)
    )
    
    # Clip to points near the tower
    # print("Clipping LAS catalog...")
    lasclip <- lidR::clip_circle(lascat, utm_coord[1, 1], utm_coord[1, 2], buffer_radius)
    
    # print("Normalizing height...")
    # Normalize using the ground points already in the cloud
    lasnorm <- lidR::normalize_height(lasclip, knnidw())
    
    #print("Saving processed LAZ...")
    lidR::writeLAS(lasnorm, file.path(outdir, paste0(site_id, ".laz")))
    
    # Delete files in the temporary directory
    #print("Deleting old files...")
    unlink(tdir, recursive=TRUE, force=TRUE) # essentially like rm -rf
  }
}
