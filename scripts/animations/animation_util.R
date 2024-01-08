library(rgl)
library(rglplus)
library(lidR)
library(sf)
library(tidyverse)
library(foreach)

# Las file must already be voxelized
make_voxel_shapelist <- function(las, resolution, foreground) {
  with <- list(x = las@data[["X"]], y = las@data[["Y"]], z = las@data[["Z"]])
  # Move cloud to the origin so there are no rgl artifacts
  with$x <- with$x - min(with$x)
  with$y <- with$y - min(with$y)
  
  n <- npoints(las)
  u <- vector("list", n)
  for(i in 1:n) {
    voxel <- rgl::cube3d(col = foreground)
    voxel <- rgl::scale3d(voxel, resolution/2, resolution/2, resolution/2)
    voxel <- rgl::translate3d(voxel , with$x[i], with$y[i], with$z[i])
    u[[i]] <- voxel
  }
  return(u)
}

# This function assumes the cloud is drawn at the origin
make_ground_bbox <- function(las, expand_factor=1.5) {
  e <- extent(las)
  x_midpoint <- (e@xmax - e@xmin) / 2
  y_midpoint <- (e@ymax - e@ymin) / 2
  
  xmin <- x_midpoint * (1 - expand_factor)
  xmax <- x_midpoint * (1 + expand_factor)
  ymin <- y_midpoint * (1 - expand_factor)
  ymax <- y_midpoint * (1 + expand_factor)
  
  # Points are ordered for drawing the quad points clockwise
  list(
    x=c(xmin, xmin, xmax, xmax),
    y=c(ymin, ymax, ymax, ymin),
    z=rep(0, 4)
  )
}

make_camera_view_vector <- function(pitch, azimuth) {
  azimuth_cartesian <- (450 - azimuth) %% 360
  return(c(
    cos(azimuth_cartesian * pi / 180), 
    sin(azimuth_cartesian * pi / 180),
    sin(pitch * pi / 180)
  ))
}

get_scene_camera_coordinate <- function(las, ref_lat, ref_lon, x_offset, 
                                        y_offset, z_offset) {
  pt <- st_sfc(st_point(c(ref_lon, ref_lat)))
  st_crs(pt) <- 4326
  pt_utm <- st_transform(pt, crs=crs(las))
  
  # Move point to origin
  pt_scene_x <- pt_utm[[1]][1] - min(las$X)
  pt_scene_y <- pt_utm[[1]][2] - min(las$Y)
  
  # Apply offset
  cam_scene_x <- pt_scene_x - x_offset
  cam_scene_y <- pt_scene_y - y_offset
  cam_scene_z <- z_offset
  
  return(c(cam_scene_x, cam_scene_y, cam_scene_z))
}

get_map_camera_coordinate <- function(las, ref_lat, ref_lon, x_offset, y_offset,
                                      z_offset) {
  pt <- st_sfc(st_point(c(ref_lon, ref_lat)))
  st_crs(pt) <- 4326
  pt_utm <- st_transform(pt, crs=crs(las))
  
  pt_utm_x <- pt_utm[[1]][1]
  pt_utm_y <- pt_utm[[1]][2]
  
  # Apply offset
  cam_utm_x <- pt_utm_x - x_offset
  cam_utm_y <- pt_utm_y - y_offset
  cam_utm_z <- z_offset
  
  return(c(cam_utm_x, cam_utm_y, cam_utm_z))
}

render_viewpoint <- function(las, vox_resolution, ref_lat, ref_lon, x_offset, y_offset,
                             z_offset, pitch, azimuth, fov, outfile, 
                             foreground=rgb(50, 168, 82, maxColorValue=255), 
                             background=rgb(77,  57, 34, maxColorValue=255),
                             save=TRUE) {
  
  not_ground <- filter_poi(las, Classification != 2)
  
  camera_map_pos <- get_map_camera_coordinate(las, ref_lat, ref_lon, x_offset,
                                              y_offset, z_offset)
  
  # Clipping points directly in front of the radiometer doesn't seem to be
  # necessary since rgl apparently applies a near clipping mask.
  
  # not_ground <- filter_poi(
  #   not_ground,
  #   (
  #     (camera_map_pos[1] - X)^2 + 
  #     (camera_map_pos[2] - Y)^2 + 
  #     (camera_map_pos[3] - Z)^2
  #   ) ^ 0.5 > vox_resolution*2
  # )
  
  vox <- voxelize_points(not_ground, vox_resolution)
  
  ground_bbox <- make_ground_bbox(vox)
  
  camera_view <- make_camera_view_vector(pitch, azimuth)
  camera_pos  <- get_scene_camera_coordinate(las, ref_lat, ref_lon,
                                             x_offset, y_offset, z_offset)
  
  shapes <- make_voxel_shapelist(vox, vox_resolution, foreground)
  
  open3d()
  shapelist3d(shapes, col=foreground)
  bg3d(color="#000000", fogtype="none")
  quads3d(ground_bbox$x, ground_bbox$y, ground_bbox$z, col=background, lit=FALSE)
  rgl.camera(position=camera_pos, direction=camera_view, fov=fov)
  par3d(windowRect=c(100, 100, 1100, 1100))
  
  if (save) {
    snapshot3d(filename=outfile, width=500, height=500, webshot=FALSE)
    close3d()
  }
}
