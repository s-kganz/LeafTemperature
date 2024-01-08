source("scripts/animations/animation_util.R")

sensor_positions <- read_csv("data_working/neon_radiometer_sensor_positions.csv") %>%
  filter(str_detect(referenceLocationID, "TOWER"),
         zOffset > 1,
         siteID %in% c("WREF", "TEAK")) 

tmax <- 10
fps <- 60
vox_resolution <- 0.25 # m
site <- "TEAK"
las <- readLAS(file.path("data_working/neon_neartower_lidar/", str_c(site, ".laz")))

idx_of_site <- first(which(sensor_positions$siteID == site))

lat <- sensor_positions$locationReferenceLatitude[idx_of_site]
lon <- sensor_positions$locationReferenceLongitude[idx_of_site]
xoffset <- sensor_positions$xOffset[idx_of_site]
yoffset <- sensor_positions$yOffset[idx_of_site]
pitch <- sensor_positions$pitch[idx_of_site]
azimuth <- sensor_positions$azimuth[idx_of_site]
zmin <- 1
zmax <- sensor_positions %>% filter(siteID == site) %>% pull(zOffset) %>% max()
zoffsets <- seq(zmin, zmax, length.out=tmax*fps)
fov <- 44 # deg, https://www.apogeeinstruments.com/content/SI-100-400-spec-sheet.pdf

view_vector <- make_camera_view_vector(pitch, azimuth)

# Animation 1: viewpoint of the radiometer moving up the tower ----
# Draw the initial scene
render_viewpoint(las, vox_resolution, lat, lon, xoffset, yoffset, zmin, pitch, 
                 azimuth, fov, "", save=FALSE)

# This function just moves the camera
frame_function <- function(t, tmax=10) {
  this_z <- approx(c(0, tmax), c(zmin, zmax), xout=t)$y
  newpos <- get_scene_camera_coordinate(las, lat, lon, xoffset, yoffset, this_z)
  
  # Position the camera
  rgl.camera(newpos, view_vector, fov=fov)
}

outdir <- "graphics"
outfile <- "teak_radiometer_perspective.mp4"
if (!file.exists(file.path(outdir, outfile)))
rgl.makemovie(
  frame_function,
  tmax=tmax,
  fps=fps,
  nframes=tmax*fps,
  output.path=outdir,
  output.filename=outfile
)

# Animation 2: wide angle view of the view cone moving up the tower ----
draw_view_cone <- function(las, lat, lon, xoffset, yoffset, zoffset, 
                           view_vector, fov=44, ...) {
  # First center of the cone
  p0 <- get_scene_camera_coordinate(las, lat, lon, xoffset, yoffset, zoffset)
  
  # Draw the cone to z = 0 or 20 m (reasonable max view distance)
  z1 <- max(0, zoffset - 20)
  
  # View vector is parameterized as z0 + view_vector*t. Determine appropriate t.
  dz <- zoffset - z1
  t <- dz / view_vector[3]
  p1 <- p0 - view_vector*t
  
  # Calculate radius at the final point to yield an overall fov of fov deg. This
  # has a half-angle of fov/2 deg.
  
  # First, height of the cone
  h <- sqrt(sum((p1 - p0)^2))
  r <- tan(fov/2 * pi / 180) * h
  
  # Now draw the cone
  cylinder3d(
    center=rbind(p0, p1),
    radius=c(0.1, r),
    ... # passed to material3d()
  )
}
cam_pos_teak <- c(-90 -50, 0) # figure out a good value manually
render_viewpoint(las, vox_resolution, lat, lon, xoffset, yoffset, zmin, pitch, 
                 azimuth, fov, "", save=FALSE)
rgl.camera(cam_pos_teak, fov=44)
