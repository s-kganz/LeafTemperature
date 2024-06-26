% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figures.R
\name{fig1_sensor_heights}
\alias{fig1_sensor_heights}
\alias{fig2_model_one_to_one}
\alias{fig3_regression_slopes}
\alias{fig4_temp_forcing}
\alias{fig5_gs_gbh_sensitivity}
\alias{fig6_shade_gpp}
\alias{fig7_le_rn_flux}
\alias{fig_s1_lidar_fit}
\alias{fig_s2_medlyn_fit}
\alias{fig_s3_aq_fit}
\alias{fig_s4_rad_model_pairings}
\alias{fig_s5_g1_sensitivity}
\alias{fig_s6_omega_density}
\alias{fig_s7_flux_boxplots}
\title{Functions to generate manuscript figures}
\usage{
fig1_sensor_heights(sensor_heights, site_meta, tower_color = "grey50")

fig2_model_one_to_one(eb_rad_tcan_summary)

fig3_regression_slopes(eb_regressions)

fig4_temp_forcing(eb_result)

fig5_gs_gbh_sensitivity(grid_eb, eb_result)

fig6_shade_gpp(shade_gpp, site_lai)

fig7_le_rn_flux(eb_result)

fig_s1_lidar_fit(lidar_constants_df, iv_io_data, max_lai = 7)

fig_s2_medlyn_fit(
  medlyn_data,
  medlyn_constants,
  gs_cutoff = 1.5,
  line_color = "grey50"
)

fig_s3_aq_fit(aq_constants, aq_data)

fig_s4_rad_model_pairings(model_radiometer_crosswalk)

fig_s5_g1_sensitivity(wref_g1_regressions)

fig_s6_omega_density(eb_result)

fig_s7_flux_boxplots(eb_result)
}
\arguments{
\item{sensor_heights}{Table of Ameriflux sensor heights as returned by \code{\link[=amf_var_heights]{amf_var_heights()}}}

\item{site_meta}{Site metadata table included with package}

\item{tower_color}{Color used for the tower pattern.}

\item{eb_rad_tcan_summary}{Table of model validation results as generated by \code{\link[=model_radiometer_comparison]{model_radiometer_comparison()}}}

\item{eb_regressions}{Table of model II regression results as generated by \code{\link[=run_neon_energy_balance]{run_neon_energy_balance()}}}

\item{eb_result}{Table of energy balance model output as generated by \code{\link[=run_neon_energy_balance]{run_neon_energy_balance()}}}

\item{grid_eb}{Table of gs/gbh sensitivity analysis as generated by \code{\link[=gs_gbh_sensitivity]{gs_gbh_sensitivity()}}}

\item{shade_gpp}{Table of proportion shade GPP in each site as generated by \code{\link[=run_neon_energy_balance]{run_neon_energy_balance()}}}

\item{lidar_constants_df, }{iv_io_data Tables generated by \code{\link[=fit_lidar_constants]{fit_lidar_constants()}}}

\item{max_lai, }{Maximum leaf area index to use in \code{\link[=fig_s1_lidar_fit]{fig_s1_lidar_fit()}}}

\item{medlyn_data, }{Medlyn_constants Tables generated by \code{\link[=fit_medlyn_slopes]{fit_medlyn_slopes()}}}

\item{gs_cutoff, }{Cutoff to exclude outliers in \code{\link[=fig_s2_medlyn_fit]{fig_s2_medlyn_fit()}}}

\item{line_color, }{Color for the 1:1 line in \code{\link[=fig_s2_medlyn_fit]{fig_s2_medlyn_fit()}}}

\item{aq_constants, }{aq_data Tables generated by \code{\link[=fit_aq_curves]{fit_aq_curves()}}}

\item{model_radiometer_crosswalk}{Crosswalk between model layer and radiometer heights generated by \code{\link[=model_radiometer_comparison]{model_radiometer_comparison()}}}

\item{wref_g1_regressions}{Model II regressions from the g1 sensitivity analysis generated by \code{\link[=run_neon_energy_balance]{run_neon_energy_balance()}}}
}
\value{
Corresponding figure from the manuscript.
}
\description{
Functions to generate manuscript figures
}
