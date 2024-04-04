library(neonUtilities)
library(lidR)
library(terra)

sites <- read.csv("data_working/neon_site_metadata.csv")
outdir <- file.path("data_working", "neon_neartower_lidar")

if (!dir.exists(outdir)) dir.create(outdir)

product <- "DP1.30003.001"
lidar_avail <- getProductInfo(product)

buffer_radius <- 50
light_sample_radius <- 10

for(i in 1:nrow(sites)) {
  site_id <- sites$site_neon[i]
  lat <- sites$tower_lat[i]
  lon <- sites$tower_lon[i]
  utm <- sites$utm_zone[i]
  
  # Find the most recent acquisition
  site_index <- match(site_id, lidar_avail$siteCodes$siteCode)
  acq_year_mo <- tail(lidar_avail$siteCodes$availableMonths[[site_index]], n=1)
  acq_year <- strsplit(acq_year_mo, "-")[[1]][1]
  
  # Convert lat lon to UTM
  utm_proj4 <- paste("+proj=utm ", "+zone=", utm, " +datum=NAD83 ", "+units=m",
                     sep="")
  utm_coord <- project(
    matrix(c(lon, lat), nrow=1), 
    from=crs("epsg:4326"), to=crs(utm_proj4)
  )
  
  utm_x <- utm_coord[1, 1]
  utm_y <- utm_coord[1, 2]
  
  # Download files to a tempdir
  tdir <- tempdir()
  
  byTileAOP(
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
  print("Reading downloaded LAS files...")
  lascat <- readLAScatalog(
    list.files(file.path(tdir, product), pattern="*.laz", recursive=TRUE,
               full.names=TRUE)
  )
  
  # Clip to points near the tower
  print("Clipping LAS catalog...")
  lasclip <- clip_circle(lascat, utm_coord[1, 1], utm_coord[1, 2], buffer_radius)
  
  print("Normalizing height...")
  # Normalize using the ground points already in the cloud
  lasnorm <- normalize_height(lasclip, knnidw())
  
  # Add a light sample attribute for calibration
  is_light_sample <- ifelse(
    (lasnorm$X - utm_x)^2 + (lasnorm$Y - utm_y)^2 < light_sample_radius^2,
    1, 0
  )
  
  lasfinal <- add_lasattribute(
    lasnorm,
    is_light_sample,
    "light_sample",
    "light_sample"
  )
  
  print("Saving processed LAZ...")
  writeLAS(lasfinal, file.path(outdir, paste0(site_id, ".laz")))
  
  # Delete files in the temporary directory
  print("Deleting old files...")
  unlink(tdir, recursive=TRUE, force=TRUE) # essentially like rm -rf
}
