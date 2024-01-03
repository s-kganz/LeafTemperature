library(lidR)
library(foreach)
library(rgl)
library(rglplus)
library(magrittr)
library(ggdark)
library(gganimate)
source("scripts/lidar_util.R")

taolist <- lapply(list(
  readTLSLAS("data_working/taos/illum_flagged/scan55_2.las"),
  readTLSLAS("data_working/taos/illum_flagged/scan42_1.las"),
  readTLSLAS("data_working/taos/illum_flagged/scan26_2.las")),
  function(las) decimate_points(las, homogenize(500))
)

plot_size <- 20
spacing   <- seq(1, 10, length.out = 20)
reps      <- 5
spacing_cross <- expand.grid(spacing=spacing, rep=1:reps)

sza   <- 30 * pi / 180 # solar zenith angle
sw_in <- 800 # incoming solar radiation

center <- c(plot_size/2, plot_size/2)
dcenter <- generate_points(2, plot_size) - c(plot_size/2, plot_size/2)

taos <- sample.int(length(taolist), nrow(dcenter), replace=TRUE)
tao_rotation <- runif(nrow(dcenter)) * 2 * pi

build_synthetic_las <- function(treepts, tao_choice, tao_rotation, tao_list,
                                header=list(), ...) {
  data <- foreach(i=1:length(tao_choice), .combine=rbind) %do% {
    this_tao <- tao_list[[tao_choice[i]]]
    this_rot <- tao_rotation[i]
    
    # rotate and offset
    this_tao <- this_tao %>% 
      rotate_las(this_rot, axis="z") %>%
      offset_las(treepts[i, 1], treepts[i, 2], xyz_only=TRUE)
  }
  LAS(data, header, ...)
}

# # Set up view parameters
# # (I am a child)
# mypp <- par3d()[c("userMatrix", "userProjection", "windowRect", "scale")]
# mypp$windowRect <- c(0, 0, 1920, 1080)
# dput(mypp, file="data_working/viewPar_synthetic.R")
pp <- dget("data_working/viewPar_synthetic.R")
pp$windowRect <- c(0, 0, 1700, 1080)
frame_func <- function(t, debug=FALSE) {
  # zoom factor is essentially equal to the frame time
  zfactor <- t/10 + 1
  newpts <- center + dcenter * zfactor
  new_las <- build_synthetic_las(newpts, taos, tao_rotation, taolist)
  
  # drop any points out of bounds
  new_las <- filter_poi(new_las, X >= 0, Y >= 0, X <= plot_size, Y <= plot_size)
  
  # rotate and get light availability
  sza <- 30 * pi / 180
  azi <- 3 * pi/2 # from right side of scene
  azi <- (7 * pi / 2 - azi) %% (2 * pi)
  
  new_rotate <- new_las %>%
    rotate_las(azi, axis="z") %>%
    rotate_las(sza, axis="y")
  
  light_avail <- light_availability(new_rotate, 0.1) * sw_in
  
  new_light <- add_attribute(new_las, light_avail, "light_availability")
  
  # lidR::plot doesn't work for custom breaks, so we have to force the color
  # scale by adding two points that serve as the min and max.
  # See https://github.com/r-lidar/lidR/issues/702
  dummy_points <- data.frame(
    X=c(20, 20),
    Y=c(20, 20),
    Z=c(0, 0),
    light_availability=c(0, sw_in)
  )
  
  new_light@data <- rbind(
    new_light@data, dummy_points, fill=TRUE
  )
  
  # Kill previous window if it is open
  if (length(rgl.dev.list())) close3d()
  # plot(wr_las, axis=TRUE, legend=TRUE)
  plot(new_light, axis=TRUE, legend=TRUE, color="light_availability",
       pal=hcl.colors(16, palette="Inferno"))
  
  text3d(20, 0, 20,
         "Shortwave radiation, W/m2",
         color="white", family="sans", cex=1.5, pos=3)
  
  # Plot a vector showing solar angle
  p1 <- c(10, 10, 15)
  len <- 5
  # draw_sun_vector(p1, len, sza, azi, color="yellow")
  par3d(pp)
  # rgl.draw()
}

tmax <- 15
fps <- 30
rgl.makemovie(frame=frame_func, tmin=0, tmax=tmax, fps=fps,
              nframes = tmax * fps,
              output.path="graphics",
              output.filename="synthetic_forest_movie.mp4",
              quiet=FALSE)

# Now animate the associated data
# rad_data <- read_csv("data_working/synthetic_forest_radiation_results.csv") %>%
#   filter(gs_forcing == 0, ta_forcing == 0)
# 
# anim <- rad_data %>%
#   mutate(
#     norm_light = light.avail.mn / sw_in,
#     z_slice = case_match(
#       z_slice,
#       0 ~ "[0.0 - 2.5] m",
#       1 ~ "[2.5 - 5.0] m",
#       2 ~ "[5.0 - 7.5] m",
#       3 ~ "> 7.5 m"
#     )
#   )  %>%
#   ggplot(aes(x=spacing, y=norm_light, color=z_slice)) + 
#   geom_point(size=3) +
#   labs(x="Average tree spacing (m)",
#        y="Relative radiation load",
#        color="Canopy position") +
#   ylim(0.6, 1) +
#   dark_theme_bw() +
#   theme(axis.text=element_text(size=24),
#         axis.title=element_text(size=32),
#         strip.text=element_text(size=24),
#         legend.text=element_text(size=24),
#         legend.title=element_text(size=32),
#         plot.title=element_text(size=32),
#         panel.grid.major = element_line(color = "grey30"),
#         panel.grid.minor = element_line(color = "grey30")) +
#   transition_time(spacing) +
#   shadow_mark()
# 
# animate(anim, width=720, height=1080, 
#         renderer=ffmpeg_renderer(format="mp4"),
#         duration=15, fps=30)
