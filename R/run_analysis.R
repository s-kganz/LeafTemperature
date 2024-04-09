run_analysis <- function(amf_user, amf_email, neon_token, work_dir, log=NULL,
                         skip="nnnnnnnnnnn") {
  # Start logging
  if (!is.null(log)) sink(log, split=TRUE)
  
  steps <- c("Download tower data", "Download flux data", "Download LiDAR",
             "Calibrate LiDAR constants", "Partition flux",
             "Interpolate meteorology", "Fit AQ curves",
             "Fit Medlyn slopes", "Run energy balance model",
             "Model vs. radiometer comparison", "Conductance sensivity analysis")
  
  skip <- (strsplit(skip, "")[[1]]) == "y"
  
  # Initial output
  cat("Output directory:", work_dir, "\n")
  cat("AMF username:", amf_user, "\n")
  cat("AMF email:", amf_email, "\n")
  cat("NEON token is", ifelse(is.null(neon_token) | is.na(neon_token), 
                              "Missing", "Provided"), "\n")
  
  # Set up directories
  if (!dir.exists(work_dir)) dir.create(work_dir)
  tower_dir <- file.path(work_dir, "amf_tower_data")
  lidar_dir <- file.path(work_dir, "neon_neartower_lidar")
  model_dir <- file.path(work_dir, "model_run")
  
  dir.create(tower_dir); dir.create(lidar_dir); dir.create(model_dir)
  
  # Load site metadata and LAI. Technically this isn't necessary but just FYI
  # these are included with the package.
  data("site_meta")
  data("manual_lai")
  
  # Download all the data we need
  if (!skip[1]) get_amf_tower_data(site_meta, work_dir, amf_user, amf_email)
  else cat("Skipping", steps[1], "\n")
  
  if (!skip[2]) {
    cat("Downloading flux takes ~90 minutes. Are you sure you want to continue?\n")
    cat("1: Continue and download flux\n")
    cat("2: Abort \n")
    choice <- readline("Your choice: ")
    if (choice == 1) get_neon_flux(site_meta, neon_token, work_dir)
    else stop("Aborted at flux download")
  } 
  else cat("Skipping", steps[2], "\n")
  
  if (!skip[3]) get_neon_lidar(site_meta, lidar_dir)
  else cat("Skipping", steps[3], "\n")
  
  # Fit lidar constants
  if (!skip[4]) fit_neon_lidar_k(site_meta, manual_lai, tower_dir, lidar_dir, work_dir)
  else cat("Skipping", steps[4], "\n")
  
  # Partition flux
  raw_flux <- read_csv(file.path(work_dir, "cross_site_flux.csv"),
                       show_col_types=FALSE)
  if (!skip[5]) partition_neon_flux(site_meta, raw_flux, tower_dir, work_dir)
  else cat("Skipping", steps[5], "\n")
  
  # Within-canopy meteorology
  partition_flux <- read_csv(file.path(work_dir, "cross_site_flux_partition_qc.csv"),
                             show_col_types=FALSE)
  lad_profiles <- read_csv(file.path(work_dir, "neon_lad_profiles.csv"),
                           show_col_types=FALSE)
  if (!skip[6]) get_within_canopy_meteorology(site_meta, manual_lai, 
                                              lad_profiles, partition_flux,
                                              work_dir)
  else cat("Skipping", steps[6], "\n")
  
  # Fit AQ curves
  if (!skip[7]) fit_aq_curves(site_meta, partition_flux, tower_dir, work_dir)
  else cat("Skipping", steps[7], "\n")
  
  
  # Fit Medlyn slopes
  if (!skip[8]) fit_medlyn_slopes(site_meta, manual_lai, partition_flux, 
                                  tower_dir, work_dir)
  else cat("Skipping", steps[8], "\n")
  
  # Run the EB model
  aq_constants <- read_csv(file.path(work_dir, "cross_site_aq_constants.csv"),
                           show_col_types=FALSE)
  lidar_constants <- read_csv(file.path(work_dir, "neon_lidar_constants.csv"),
                              show_col_types=FALSE)
  medlyn_constants <- read_csv(file.path(work_dir, "cross_site_medlyn_coefficients.csv"),
                               show_col_types=FALSE)
  interp_meteo <- read_csv(file.path(work_dir, "cross_site_interpolated_meteorology.csv"),
                           show_col_types=FALSE)
  
  if (!skip[9])
    run_neon_energy_balance(
      site_meta,
      manual_lai,
      partition_flux,
      aq_constants,
      medlyn_constants,
      lidar_constants,
      interp_meteo,
      tower_dir,
      model_dir
    )
  else cat("Skipping", steps[9], "\n")
  
  eb_result <- read_csv(file.path(model_dir, "cross_site_eb.csv"),
                        show_col_types=FALSE)
  rad_tcan <- read_csv(file.path(work_dir, "cross_site_tcan.csv"),
                       show_col_types=FALSE)
  if (!skip[10]) model_radiometer_comparison(eb_result, rad_tcan, model_dir)
  else cat("Skipping", steps[10], "\n")
  
  if (!skip[11]) gs_gbh_sensitivity(model_dir)
  else cat("Skipping", steps[11], "\n")
  
  # Reset logging
  if (!is.null(log)) sink()
}
