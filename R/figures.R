# This file provides functions to recreate figures in the paper. By default,
# all the necessary data gets loaded instead of just the data needed
# for a particular figure.

library(tidyverse)
library(ggpattern)
library(amerifluxr)
library(photosynthesis)
library(lmodel2)
library(metR)

# Main figures ----
fig1_sensor_heights <- function(sensor_heights, site_meta, tower_color="grey50") {
  # Filter sensor heights to radiometers and within-canopy microclimate
  sensor_regex <- paste0(c("T_CANOPY_\\d_\\d_\\d", "H2O_\\d_\\d_\\d",
                           "TA_\\d_\\d_\\d", "WS_\\d_\\d_\\d"),
                         collapse="|")
  
  sensor_heights_filter <- sensor_heights %>%
    filter(str_detect(Variable, sensor_regex)) %>%
    select(Site_ID, Variable, Height) %>%
    mutate(site_neon = site_meta$site_neon[match(Site_ID, site_meta$site_ameriflux)],
           Height = floor(Height/2)*2,
           vartype = ifelse(str_detect(Variable, "T_CANOPY"), "Radiometer", "Microclimate\nTriplet")) %>%
    select(-Variable) %>%
    distinct() %>%
    # Drop radiometers below the canopy
    filter(vartype != "Radiometer" | Height > 2)
  
  # Attach filepath of the tree svgs to metadata table
  site_meta_filepath <- site_meta %>%
      mutate(filepath = file.path(
      "graphics", "tree_svgs", str_c(str_to_lower(site_ameriflux), ".png")
    ))
  
  # Make a named vector to match pngs with site
  file_match_vector <- site_meta_filepath$filepath
  names(file_match_vector) <- site_meta_filepath$site_ameriflux
  
  ggplot(site_meta_filepath, aes(x=site_ameriflux, y=tower_height)) +
    # Add the tower top first, this can make it easier to not clip
    # the tower top if the image is resized.
    geom_bar_pattern(
      aes(x="right", y=tower_height),
      fill="black",
      width=0.3,
      stat="identity",
      pattern="image",
      pattern_type = "none",
      pattern_filename="graphics/tower_top.png",
      pattern_scale=-1,
      pattern_gravity="north",
      pattern_yoffset=0.5,
      alpha=0
    ) +
    # Tower tile
    geom_bar_pattern(
      aes(x="right", y=tower_height*0.8),
      stat="identity",
      fill="black",
      pattern="image",
      pattern_type="tile",
      pattern_filename="graphics/tower_tile.png",
      pattern_filter="box",
      pattern_scale=-1,
      color=tower_color,
      linewidth=1,
      width=0.2,
      alpha=0
    ) +
    # Now add the trees
    geom_bar_pattern(
      aes(pattern_filename=site_ameriflux, y=canopy_height, x="left"),
      stat="identity",
      pattern="image",
      pattern_gravity="south",
      alpha=0
    ) +
    scale_pattern_filename_manual(values=file_match_vector) +
    # Debug geom to make sure trees match the actual height
    # geom_point(aes(y=canopy_height, x="left"), color="green") +
    geom_point(
      data=sensor_heights_filter,
      mapping=aes(x="right", y=Height, fill=vartype, shape=vartype),
      position=position_dodge(width=0.5),
      col="black",
      size=2.5
    ) +
    scale_shape_manual(values=c(21, 25)) +
    geom_label(aes(label=site_neon), x=-Inf, y=Inf, vjust=1, hjust=0) +
    facet_wrap(~ site_neon, scales="free_y") +
    # suppress image legend
    guides(pattern_filename="none") +
    # Styling
    labs(x="", y="Height above ground (m)", fill="", shape="") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_blank(),
          legend.text = element_text(size=12),
          panel.spacing = unit(1.5, "lines"),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14))
}
fig2_model_one_to_one <- function(eb_rad_tcan_summary) {
  p <- ggplot(eb_rad_tcan_summary, aes(x=TCAN, y=EB_MODEL_Tl)) +
    geom_bin2d() +
    geom_abline(slope=1, intercept=0, color="black", linetype="dashed") +
    scale_fill_viridis_c(limits=c(0, 200), oob=scales::squish) +
    facet_wrap(~ SITE_NEON) +
    theme_bw() +
    coord_equal() +
    labs(x=expression("Radiometer temperature ("*degree*"C)"),
         y=expression("Model temperature ("*degree*"C)"),
         fill="Count") +
    theme(legend.direction = "horizontal")

  lemon::reposition_legend(p, "center", panel=c("panel-3-2", "panel-3-3"))
}
fig3_regression_slopes <- function(eb_regressions) {
  eb_regressions %>%
    mutate(LAYER_L = factor(LAYER_L)) %>%
    ggplot(aes(x=slope_p50, y=LAYER_L)) +
    geom_vline(xintercept=1, color="grey50", linetype="dashed") +
    geom_errorbar(aes(xmin=slope_p025, xmax=slope_p975)) +
    geom_point(aes(fill=prop_Tl_lt_Ta), color="black", size=4, pch=21) +
    geom_label(aes(label=SITE_NEON), x=-Inf, y=Inf,
               hjust="left", vjust="top", size=5) +
    scale_fill_viridis_c() +
    scale_y_discrete(limits=rev) +
    scale_x_continuous(limits=c(0.97, 1.15)) +
    facet_grid(SITE_NEON ~ .,
               scales="free_y", space="free") +
    theme_bw() +
    theme(strip.placement="outside",
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          axis.title = element_text(size=18),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14)
    ) +
    guides(fill = guide_colorbar(
      ticks.colour="black",
      frame.colour="black"
    )) +
    labs(x=expression("T"[L]~"vs"~"T"[A]~"regression slope ("~degree*C*"/"*degree*C*")"),
         y=expression("Cumulative LAI (m"^2~"m"^-2*")"),
         fill=expression("Proportion T"[L]~"<"~"T"[A]))
}
fig4_temp_forcing <- function(eb_result) {
  eb_result_summary <- eb_result %>%
    group_by(SITE_NEON) %>%
    mutate(
      canopy_group = case_when(
        LAYER_L == max(LAYER_L) ~ "Lower",
        LAYER_L == 0 ~ "Upper",
        .default="Middle"
      ),
      canopy_group = factor(canopy_group, levels=c("Lower", "Middle", "Upper"))
    ) %>%
    group_by(SITE_NEON, canopy_group) %>%
    summarize(
      Rn_forcing = mean(EB_MODEL_dT_Rn),
      Rn_sd = sd(EB_MODEL_dT_Rn),
      LE_forcing = -mean(EB_MODEL_dT_LE),
      LE_sd = sd(EB_MODEL_dT_LE),
      net = Rn_forcing + LE_forcing
    )
  
  p <- ggplot(eb_result_summary, aes(y=canopy_group)) +
    geom_bar(aes(x=Rn_forcing, fill="Rn"), stat="identity") +
    geom_bar(aes(x=LE_forcing, fill="LE"), stat="identity") +
    geom_errorbar(aes(x=Rn_forcing, xmin=Rn_forcing-Rn_sd, 
                      xmax=Rn_forcing+Rn_sd),
                  width=0.3) +
    geom_errorbar(aes(x=LE_forcing, xmin=LE_forcing-LE_sd, 
                      xmax=LE_forcing+LE_sd),
                  width=0.3) +
    # geom_point(aes(x=net)) +
    facet_wrap(~ SITE_NEON) +
    scale_fill_manual(
      labels=c("Transpiration", "Net radiation"),
      values=hcl(h=c(195, 15), l=65, c=100)
    ) +
    theme_bw() +
    labs(x="Temperature forcing (K)",
         y="", fill="")
  
  lemon::reposition_legend(p, 'center', panel=c("panel-3-3"))
}
fig5_gs_gbh_sensitivity <- function(grid_eb) {
  ggplot(grid_eb) +
    geom_contour_fill(mapping=aes(x=gs, y=gbH, z=dT),
                      binwidth=0.1) +
    # geom_contour(mapping=aes(x=gs, y=gbH, z=dT), color="black",
    #              binwidth=0.25) +
    scale_fill_divergent() +
    # geom_label_contour(mapping=aes(x=gs, y=gbH, z=dT),
    #                    skip=0, label.placer=label_placer_fraction()) +
    facet_wrap(~ tot_rad,
               labeller=label_bquote(LW[abs]*"+"*SW[abs]~"="~.(tot_rad)~"W m"^-2)) +
    #xlim(c(gs_min, gs_max)) + ylim(gbH_min, gbH_max) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    theme_bw() +
    labs(x=expression("Stomatal conductance"~"(mol m"^-2~"s"^-1*")"),
         y=expression("Boundary layer conductance"~"(mol m"^-2~"s"^-1*")"),
         fill=expression(Delta*"T (K)"))
}
fig6_shade_gpp <- function(shade_gpp) {
  # Shaded GPP calculation
  shade_gpp_sorted <- shade_gpp %>% arrange(desc(prop_gpp_shade)) %>% 
    pull(SITE_NEON)
  
  shade_gpp %>%
    ggplot(aes(y=SITE_NEON, x=prop_gpp_shade)) +
    geom_point() +
    geom_vline(xintercept=0.552, color="red", linetype="dashed") +
    scale_y_discrete(limits=shade_gpp_sorted) +
    theme_bw() +
    labs(x="Proportion shaded GPP", y="") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank()) +
    xlim(0, 0.7)
}


# Supplemental figures ----
fig_s1_lidar_fit <- function(lidar_constants_df, iv_io_data, 
                             max_lai=7) {
  clai <- seq(0, max_lai, length.out=20)
  iv_io_extrapolate <- expand.grid(clai, 1:nrow(lidar_constants_df)) %>%
    cbind(lidar_constants_df[.$Var2, ]) %>%
    select(-Var2) %>% rename(clai = Var1) %>%
    mutate(Beam = exp(kb_real * clai),
           Diffuse = exp(kd_real * clai)) %>%
    pivot_longer(c(Beam, Diffuse))
  
  ggplot(iv_io_data, aes(x=cum_lai, y=iv_io_mean)) +
    geom_point(aes(color=proportion_diffuse_bin)) +
    geom_line(data=iv_io_extrapolate, 
              mapping=aes(x=clai, y=value, group=name, linetype=name)) +
    facet_wrap(~ site, scales="free") +
    labs(x="Cumulative LAI", y="Proportion available light",
         color="Proportion diffuse light", linetype="") +
    scale_y_log10() + theme_bw()
}
fig_s2_medlyn_fit <- function(medlyn_data, medlyn_constants, 
                              gs_cutoff=1.5, line_color="grey50") {
  # Some of the P-M Gs values are unreasonable. This are included in fitting
  # but should be excluded from the figure for clarity.
  medlyn_data %>%
    ggplot(aes(x=Gs_mol, y=Gs_mol_predicted)) +
    geom_point(aes(color=VPD)) +
    geom_label(
      data=medlyn_constants,
      mapping=aes(x=0, y=gs_cutoff, label=str_c("g1 = ", format(round(g1, 2), nsmall=2))),
      vjust=1,
      hjust=0
    ) +
    geom_abline(slope=1, intercept=0, color=line_color) +
    coord_equal(xlim=c(0, gs_cutoff), ylim=c(0, gs_cutoff)) +
    scale_color_viridis_c() +
    facet_wrap(~ site_neon) +
    theme_bw() +
    labs(x=expression("Penman-Monteith G"[s]~"(mol m"^-2~"s"^-1*")"),
         y=expression("Model G"[s]~"(mol m"^-2~"s"^-1*")"),
         color="VPD (kPa)")
}
fig_s3_aq_fit <- function(aq_constants, aq_data) {
  # Show fit on top of raw data
  sw_range <- seq(0, 1250)
  fit_pred <- aq_constants %>%
    mutate(
      A_pred = pmap(
        list(k_sat, phi_J, theta_J, Rd),
        function(k_sat, phi_J, theta_J, Rd) {
          data.frame(
            A = marshall_biscoe_1980(sw_range, k_sat, phi_J, theta_J) - Rd,
            Q = sw_range
          )
        }
      )
    ) %>%
    select(site, A_pred) %>% unnest(A_pred)
  
  ggplot(NULL) +
    geom_point(mapping=aes(x=.Qabs, y=.A), data=aq_data, alpha=0.1) +
    geom_line(mapping=aes(x=Q, y=A, color="Fit"), data=fit_pred, lwd=1, show.legend = FALSE) +
    facet_wrap(~ site, scales="free_y") +
    labs(x=expression("Top-of-canopy shortwave (W m"^2*")"),
         y=expression("Canopy GPP ("*mu*"mol CO"[2]~"m"^-2~s^-1*")"),
         color="") +
    theme_bw()
}
fig_s4_rad_model_pairings <- function(model_radiometer_crosswalk) {
  # This approach actually duplicates some points on the radiometer side. But
  # you can't see this in the final product so whatever.
  ggplot(model_radiometer_crosswalk) +
    geom_segment(aes(x="Model Layers", xend="Radiometers", y=LAYER_z, yend=Height),
                 linetype="dashed") +
    geom_point(aes(x="Model Layers", y=LAYER_z), fill="orange", pch=22, size=2) +
    # Match color and pch in fig 1
    geom_point(aes(x="Radiometers", y=Height), fill="#00BFC4", pch=25, size=2) +
    facet_wrap(~ SITE_NEON, scales="free_y") +
    theme_bw() +
    labs(x="", y="Height above ground (m)")
}
fig_s5_g1_sensitivity <- function(wref_g1_regressions) {
  ggplot(wref_g1_regressions, aes(y=factor(LAYER_L))) +
    geom_vline(xintercept=1, color="grey50", linetype="dashed") +
    geom_errorbar(aes(xmin=slope_p025, xmax=slope_p975)) +
    geom_point(aes(x=slope_p50, fill=prop_Tl_lt_Ta), size=2, pch=21) +
    facet_wrap(~ str_c("g1 = ", MEDLYN_g1)) +
    scale_fill_viridis_c() +
    scale_y_discrete(limits=rev) +
    scale_x_continuous(limits=c(0.97, NA)) +
    theme_bw() +
    labs(x=expression("T"[leaf]~"vs"~"T"[air]~"regression slope ("~degree*C*"/"*degree*C*")"),
         y=expression("Cumulative LAI (m"^2~"m"^-2*")"),
         fill=expression("Proportion T"[leaf]~"<"~"T"[air])) +
    guides(fill = guide_colorbar(
      ticks.colour="black",
      frame.colour="black"
    ))
}

# Spit out all figures ----
# Helper function to prevent overwrites
safe_save <- function(filename, plot, allow_overwrite=FALSE, ...) {
  if (allow_overwrite | !file.exists(filename)) {
    ggsave(filename, plot, ...)
  } else {
    warning(filename, " already exists! Figure not written.")
  }
}

write_all_figures <- function(site_meta, search_dir, out_dir, overwrite=FALSE) {
  # Read data ----
  use_sites <- c("ABBY", "DEJU", "JERC", "OSBS", "RMNP", "TALL", "WREF")
  
  # General site metadata
  site_meta <- site_meta %>%
    filter(site_neon %in% use_sites)
  
  # Height of ameriflux sensors
  sensor_heights <- amf_var_heights() %>%
    filter(Site_ID %in% site_meta$site_ameriflux) %>%
    mutate(SITE_NEON = site_meta$site_neon[match(Site_ID, site_meta$site_ameriflux)])
  
  # Lidar constants and mean transmission data
  lidar_constants_df <- read_csv(file.path(search_dir, "neon_lidar_constants.csv")) %>%
    filter(site %in% use_sites)
  
  iv_io_data <- read_csv(file.path(search_dir, "neon_lidar_transmission_data.csv")) %>%
    filter(site %in% use_sites)
  
  # Medlyn model coefficients and fit
  medlyn_data <- read_csv(file.path(search_dir, "cross_site_medlyn_fit_results.csv")) %>%
    mutate(site_neon = site_meta$site_neon[match(site, site_meta$site_ameriflux)])
  
  medlyn_constants <- read_csv(file.path(search_dir, "cross_site_medlyn_coefficients.csv")) %>%
    mutate(site_neon = site_meta$site_neon[match(site, site_meta$site_ameriflux)]) %>%
    drop_na()
  
  # AQ curve coefficients and fit
  aq_constants <- read_csv(file.path(search_dir, "cross_site_aq_constants.csv")) %>%
    filter(site %in% use_sites)
  
  aq_data <- read_csv(file.path(search_dir, "cross_site_aq_data.csv")) %>%
    filter(site %in% use_sites)
  
  # Model layer and radiometer pairings
  model_radiometer_crosswalk <- read_csv(file.path(search_dir, "model_run", "model_rad_height_crosswalk.csv")) %>%
    # Attach NEON site code for consistency
    mutate(SITE_NEON = site_meta$site_neon[match(SITE_AMF, site_meta$site_ameriflux)])
  
  # EB model result and lmodel2 regressions
  eb_result <- read_csv(file.path(search_dir, "model_run", "cross_site_eb.csv"))
  eb_regressions <- read_csv(file.path(search_dir, "model_run", "cross_site_tl_ta_slopes.csv"))
  
  # Layer-by-layer model and radiometer comparison
  eb_rad_tcan_summary <- read_csv(file.path(search_dir, "model_run", "cross_site_model_rad_comparison.csv"))
  eb_rad_height_xwalk <- read_csv(file.path(search_dir, "model_run", "model_rad_height_crosswalk.csv"))
  
  # gs/gbh sensitivity analysis
  grid_eb <- read_csv(file.path(search_dir, "model_run", "gs_gbh_sensitivity.csv"))
  
  # Shaded GPP calculation
  shade_gpp <- read_csv(file.path(search_dir, "model_run", "cross_site_shade_gpp.csv"))
  
  # g1 sensitivity analysis
  wref_g1_sensitivity <- read_csv(file.path(search_dir, "model_run", "wref_eb_g1_sensitivity.csv"))
  wref_g1_regressions <- read_csv(file.path(search_dir, "model_run", "wref_g1_tl_ta_slopes.csv"))
  
  ## Generate all the figs ----
  ### Main text ----
  safe_save(file.path(out_dir, "fig1_sensor_heights.png"),
            fig1_sensor_heights(sensor_heights, site_meta),
            allow_overwrite=overwrite,
            width=6.5, height=4)
  
  safe_save(file.path(out_dir, "fig2_model_one_to_one.png"),
            fig2_model_one_to_one(eb_rad_tcan_summary),
            allow_overwrite=overwrite,
            width=5, height=4.5)
  
  safe_save(file.path(out_dir, "fig3_regression_slopes.png"),
            fig3_regression_slopes(eb_regressions),
            allow_overwrite=overwrite,
            width=8, height=10)
  
  safe_save(file.path(out_dir, "fig4_temp_forcing.png"),
            fig4_temp_forcing(eb_result),
            allow_overwrite=overwrite,
            width=5, height=4.5)
  
  safe_save(file.path(out_dir, "fig5_gs_gbh_sensitivity.png"),
            fig5_gs_gbh_sensitivity(grid_eb),
            allow_overwrite=overwrite,
            width=6.5, height=3)
  
  safe_save(file.path(out_dir, "fig6_shade_gpp.png"),
            fig6_shade_gpp(shade_gpp),
            allow_overwrite=overwrite,
            width=4, height=2)
  
  ### Supplemental ----
  safe_save(file.path(out_dir, "fig_s1_lidar_fit.png"),
            fig_s1_lidar_fit(lidar_constants_df, iv_io_data),
            allow_overwrite=overwrite,
            width=8, height=6)
  
  safe_save(file.path(out_dir, "fig_s2_medlyn_fit.png"),
            fig_s2_medlyn_fit(medlyn_data, medlyn_constants),
            allow_overwrite=overwrite,
            width=6, height=6)
  
  safe_save(file.path(out_dir, "fig_s3_aq_fit.png"),
            fig_s3_aq_fit(aq_constants, aq_data),
            allow_overwrite=overwrite,
            width=8, height=6)
  
  safe_save(file.path(out_dir, "fig_s4_rad_model_pairings.png"),
            fig_s4_rad_model_pairings(model_radiometer_crosswalk),
            allow_overwrite=overwrite,
            width=8, height=6)
  
  safe_save(file.path(out_dir, "fig_s5_g1_sensitivity.png"),
            fig_s5_g1_sensitivity(wref_g1_regressions),
            allow_overwrite=overwrite,
            width=8, height=6)
}
