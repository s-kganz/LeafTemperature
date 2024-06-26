% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_neon_lidar_k.R
\name{fit_neon_lidar_k}
\alias{fit_neon_lidar_k}
\title{Fit empirical light attenuation coefficients to aerial LiDAR}
\usage{
fit_neon_lidar_k(site_meta, site_lai, tower_dir, lidar_dir, outdir)
}
\arguments{
\item{site_meta}{Site metadata table \code{\link{site_meta}}.}

\item{site_lai}{Site leaf area index table \code{\link{manual_lai}}}

\item{tower_dir}{Directory containing tower data files generated by \code{\link[=get_amf_tower_data]{get_amf_tower_data()}}.}

\item{lidar_dir}{Directory containing aerial LiDAR clouds generated by \code{\link[=get_neon_lidar]{get_neon_lidar()}}}

\item{outdir}{Directory where output tables will be written.}
}
\description{
To run the radiation submodel, we need empirical light attenutation coefficients for direct and diffuse radiation. These are calculated empirically based on observations where diffuse and direct light dominate. We also need to map vertical heights above ground to cumulative leaf area index. This script calculates the attenuation coefficients and generates LAI-Z crosswalk tables.
}
