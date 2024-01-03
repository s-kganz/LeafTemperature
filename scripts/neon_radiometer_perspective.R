source("scripts/animation_util.R")

# mylas <- readLAS("data_working/neon_neartower_lidar/TEAK.laz")
# ref_lat <-   37.005842
# ref_lon <- -119.006039
# x_offset <- 2.41
# y_offset <- 5.88
# z_offset <- 10.11
# pitch <- -68
# azimuth <- 180
# fov <- 44
# 
# render_viewpoint(las, 0.25, ref_lat, ref_lon, x_offset, y_offset, z_offset,
#                  pitch, azimuth, fov, outfile="test_render.png")

out_dir <- file.path("data_working", "neon_radiometer_perspectives")

sensor_positions <- read_csv("data_working/neon_radiometer_sensor_positions.csv") %>%
  filter(str_detect(referenceLocationID, "TOWER"),
         zOffset > 1)

for(i in 1:nrow(sensor_positions)) {
  print(i)
  
  rowlist <- as.list(sensor_positions[i, ])
  out_f_name <- str_c(rowlist$siteID, rowlist$zOffset, ".png", sep="_")
  
  this_las <- readLAS(
    file.path("data_working", "neon_neartower_lidar", str_c(rowlist$siteID, ".laz"))
  )
  
  render_viewpoint(
    this_las, 0.25, rowlist$locationReferenceLatitude, rowlist$locationReferenceLongitude,
    rowlist$xOffset, rowlist$yOffset, rowlist$zOffset, rowlist$pitch,
    rowlist$azimuth, 44, outfile=file.path(out_dir, out_f_name)
  )
  
  print(file.path(out_dir, out_f_name))
}
