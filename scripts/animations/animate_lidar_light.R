# This file makes an animation showing light availability throughout the WREF
# canopy for a representative day in August, 2021.

library(lidR)
library(rgl)
library(rglplus)
library(suncalc)
library(hms)
source("scripts/lidar_util.R")
#source("scripts/download_wref_flux.R")

# Load source data
wr_las <- readLAS("data_working/neon_neartower_lidar/WREF.laz")

cloud_center_x <- mean(wr_las$X)
cloud_center_y <- mean(wr_las$Y)

wr_flux <- read_csv(file.path("data_in", "wr_flux", "wr_flux.csv")) %>%
  filter(year(TIMESTAMP) == 2021, month(TIMESTAMP) == 8,
         hour(TIMESTAMP) >= 6, hour(TIMESTAMP) <= 21) %>%
  rowwise() %>%
  mutate(SW_IN = mean(c_across(matches("SW_IN_\\d_\\d_\\d")), na.rm=TRUE)) %>%
  select(TIMESTAMP, SW_IN) %>%
  mutate(hms = as_hms(TIMESTAMP)) %>%
  group_by(hms) %>%
  summarize(SW_IN = pmax(0, median(SW_IN, na.rm=TRUE))) %>%
  mutate(dt = as.Date("2021-08-01") + as.duration(hms)) %>%
  filter(SW_IN > 0)

wr_lat <- 45.82049
wr_lon <- -121.95191

# Tweening functions for between measurements
dt_approx <- approxfun(1:nrow(wr_flux), wr_flux$dt, rule=2)
sw_approx <- approxfun(1:nrow(wr_flux), wr_flux$SW_IN, rule=2)

# Animation parameters
anim_time <- 15 # sec
fps <- 30
pp <- dget("data_working/viewPar_light.R")
pp$windowRect <- c(0, 0, 1080, 720)
pp$zoom <- 0.9

frame_func <- function(t, debug=FALSE) {
  # Calculate solar angle
  dt_t <- dt_approx(t / anim_time * nrow(wr_flux))
  sw_t <- max(1e-3, sw_approx(t / anim_time * nrow(wr_flux)))
  
  dt_t <- as_datetime(dt_t) %>% force_tz("US/Pacific")
  
  sunpos <- getSunlightPosition(date=dt_t, lat=wr_lat, lon=wr_lon)
  
  zenith <- pi / 2 - sunpos$altitude
  azimuth <- sunpos$azimuth
  # Convert azimuth to cartesian angle
  azimuth <- (7 * pi / 2 - azimuth) %% (2 * pi)
  
  if (debug) { 
    cat(strftime(dt_t), "\n")
    cat(zenith * 180 / pi, "\n")
    cat(azimuth * 180 / pi, "\n")
    cat(sw_t, "\n")
  }
  
  # Get light availability at this angle - handle edge case where sun is below
  # the horizon as all zeros
  if (zenith > (pi / 2)) {
    light_avail <- rep(NA, nrow(wr_las))
  } else {
    wr_rotate <- wr_las %>%
      rotate_las(azimuth, axis="z") %>%
      rotate_las(zenith, axis="y")
    
    light_avail <- light_availability(wr_rotate) * sw_t
  }

  # Add attribute
  wr_light <- add_attribute(wr_las, light_avail, "light_availability")
  
  # lidR::plot doesn't work for custom breaks, so we have to force the color
  # scale by adding two points that serve as the min and max.
  # See https://github.com/r-lidar/lidR/issues/702
  min_light <- 0
  max_light <- 900
  dummy_points <- data.frame(
    X=rep(cloud_center_x, 2),
    Y=rep(cloud_center_y, 2),
    Z=c(0, 0),
    light_availability=c(min_light, max_light)
  )
  
  wr_light@data <- rbind(
    wr_light@data, dummy_points, fill=TRUE
  )
  
  # Plot the points
  # rgl.hold()
  close3d()
  # plot(wr_las, axis=TRUE, legend=TRUE)
  plot(wr_light, axis=TRUE, legend=TRUE, color="light_availability",
       pal=hcl.colors(16, palette="Inferno"))
  
  text3d(160, 0, 40,
         strftime(dt_t), 
         color="white", family="sans", cex=1.5, pos=3)
  
  text3d(160, 0, 30,
         "Shortwave radiation, W/m2", 
         color="white", family="sans", cex=1.5, pos=3)

  # Plot a vector showing solar angle
  p1 <- c(50, 50, 80)
  len <- 25
  draw_sun_vector(p1, len, zenith, azimuth, color="yellow")
  par3d(pp)
  # rgl.draw()
}

rgl.makemovie(frame=frame_func, tmin=0, tmax=anim_time, fps=fps,
              nframes = anim_time * fps,
              output.path="graphics",
              output.filename="test_render.mp4",
              quiet=FALSE)
