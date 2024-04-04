library(lidR)
library(tidyverse)

rotation_matrix_x <- function(theta) {
  # Generate a rotation matrix of theta radians in 3D about the x-axis
  # (points move CCW in the Y-Z plane)
  matrix(
    c(1, 0, 0,
      0, cos(theta), -sin(theta),
      0, sin(theta), cos(theta)),
    nrow=3, ncol=3
  )
}

rotation_matrix_y <- function(theta) {
  # Generate a rotation matrix of theta radians in 3D about the x-axis
  # (points move CCW in the X-Z plane)
  matrix(
    c(
       cos(theta), 0, sin(theta),
       0,          1, 0,
      -sin(theta), 0, cos(theta)
    ),
    nrow=3, ncol=3
  )
}

rotation_matrix_z <- function(theta) {
  # Generate a rotation matrix of theta radians in 3D about the z-axis
  # (points move CW in the X-Y plane).
  matrix(
    c(cos(theta), -sin(theta), 0,
      sin(theta),  cos(theta), 0,
               0,           0, 1),
    nrow=3, ncol=3
  )
}

rotate_points <- function(data, rotation_matrix) {
  t(rotation_matrix %*% t(data))
}

rotate_las <- function(las, theta, axis="x") {
  if (axis == "x") {
    rotation_matrix <- rotation_matrix_x(theta)
  } else if (axis == "y") {
    rotation_matrix <- rotation_matrix_y(theta)
  } else {
    rotation_matrix <- rotation_matrix_z(theta)
  }
  
  las_data <- as.matrix(las@data[, 1:3])
  las_data_rotated <- rotate_points(las_data, rotation_matrix)
  
  las_copy <- las
  las_copy@data[, 1] <- las_data_rotated[, 1]
  las_copy@data[, 2] <- las_data_rotated[, 2]
  las_copy@data[, 3] <- las_data_rotated[, 3]
  
  return(las_copy)
}

distance_threshold <- function(las, pct, origin=c(0, 0)) {
  # Given a LAS file and a percentage, find the minimum distance from the origin
  # of the file where pct% of points are present.
  d <- ((las$X-origin[1])^2 + (las$Y-origin[2])^2) ^ 0.5
  this_cdf <- ecdf(d)
  sample_pts <- seq(min(d), max(d), length.out=100)
  return(
    sample_pts[which.min(abs(this_cdf(sample_pts) - pct))]
  )
}

stand_level_rumple <- function(las, res=0.5) {
  canopy_top <- decimate_points(las, highest(res=res))
  return(
    rumple_index(canopy_top$X, canopy_top$Y, canopy_top$Z)
  )
}

stand_pct_cover <- function(las, res=0.05) {
  las_count  <- grid_metrics(las, ~length(X) > 0, res=res)
  plot_rad   <- (extent(las_count)@xmax - extent(las_count)@xmin) / 2
  plot_rad_px <- plot_rad / res
  cover <- sum(las_count$V1@data@values, na.rm=TRUE) / (pi * plot_rad_px^2)
  return(cover)
}

height_threshold <- function(las, search_interval=c(0.5, 10)) {
  # Finds the minimum height in the point cloud separating understory and
  # overstory vegetation. Search the first derivative of the probability 
  # density function for heights in the point cloud for a zero in the search
  # interval.
  las_density <- density(las$Z)
  d_savgol  <- pracma::savgol(las_density$y, 7, dorder=1)
  d_roots   <- rootSolve::uniroot.all(
    approxfun(las_density$x, d_savgol), interval=search_interval
  )
  if (length(d_roots) == 0) {
    warning("[height_threshold] No threshold found, returning lower bound
            of interval")
    return(search_interval[1])
  } else {
    return(d_roots[1])
  }
}

normalize_las <- function(las, cls_algo=pmf(1, 0.1), raster_res=1,
                          nrm_algo = knnidw()) {
  ground <- classify_ground(las, cls_algo)
  dtm  <- rasterize_terrain(ground, res = raster_res, algorithm = nrm_algo)
  return(normalize_height(ground, dtm = dtm, algorithm = nrm_algo))
}

light_availability <- function(las, quadrat_res=1, kI=0.0375) {
  las@data %>% as.data.frame() %>%
    mutate(
      quad_x = round(X/quadrat_res) * quadrat_res,
      quad_y = round(Y/quadrat_res) * quadrat_res
    ) %>%
    group_by(quad_x, quad_y) %>%
    mutate(
      z_rank = rank(Z),
      dz     = max(Z) - Z,
      light  = exp(-kI * dz * log(n() / z_rank))
    ) %>%
    pull(light)
}

light_availability2 <- function(las, quadrat_res=1, k=0.0375) {
  las@data %>% as.data.frame() %>%
    mutate(
      quad_x = round(X/quadrat_res) * quadrat_res,
      quad_y = round(Y/quadrat_res) * quadrat_res
    ) %>%
    group_by(quad_x, quad_y) %>%
    mutate(
      z_rank = rank(Z),
      light  = (z_rank / n()) ^ k
    ) %>%
    pull(light)
}

sample_light_availability <- function(las, light_avail, z, slice_size=0.5) {
  # Get geometric mean light availability at each of the provided z slices
  sapply(z, function(z) {
    zmin <- z - slice_size
    zmax <- z + slice_size
    lidar_subset <- las$Z > zmin & las$Z < zmax
    
    if (!any(lidar_subset)) {
      warning("[sample_light_availability] No points found in Z range! Returning NA")
      return(NA)
    }
    
    # Geometric mean light availability in each subset
    light_avail[lidar_subset] %>% log() %>% mean() %>% exp()
  })
}

sample_attribute <- function(las, attr, subset, stat_fn) {
  stat_fn(las[[attr]][subset])
}

draw_sun_vector <- function(p1, len, zenith, azimuth, ...) {
  p0 <- c(
    p1[1] + len * cos(azimuth),
    p1[2] + len * sin(azimuth),
    p1[3] + len * cos(zenith))
  rgl::arrow3d(p0, p1, ...)
}

#' Generate 2D points in an area with given spacing.
#' 
#' Equivalent to generating points with `density = spacing^-2`.
#'
#' @param spacing Minimum distance between points
#' @param size Size of the square in which to generate points
#'
#' @return
#' @export
#'
#' @examples
generate_points <- function(spacing, size) {
  # Initialize points
  pts_to_generate <- round(10 * spacing^-2 * size ^ 2)
  X <- matrix(runif(pts_to_generate * 2), ncol=2) * size
  
  rejected_pts <- c()
  accepted_pts <- c()
  candidates   <- 1:nrow(X)
  
  while (length(candidates) > 0) {
    next_pt <- sample(candidates, 1)
    
    accepted_pts <- c(accepted_pts, next_pt)
    
    dist_from_accepted <- sqrt(apply(sweep(X, 2, X[next_pt, ])^2, 1, sum))
    
    to_reject <- which(dist_from_accepted < spacing)
    
    candidates <- setdiff(candidates, to_reject)
  }
  
  return(X[accepted_pts, ])
}

#' Add a constant x/y/z offset to a LAS dataset
#'
#' @param las LAS dataset
#' @param dx offset in x direction
#' @param dy offset in y direction
#' @param dz offset in z direction
#' @param xyz_only only return the xyz columns of the data matrix?
#'
#' @return
#' @export
#'
#' @examples
offset_las <- function(las, dx, dy, dz=0, xyz_only=FALSE) {
  las_copy <- las
  las_copy@data[, 1:3] <- sweep(las@data[, 1:3], 2, c(dx, dy, dz), "+")
  
  if (xyz_only) {
    return(las_copy@data[, 1:3])
  } else {
    return(las_copy)
  }
}
