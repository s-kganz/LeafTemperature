run_analysis <- function(amf_user, amf_email, neon_token, work_dir, log=NULL,
                         skip=rep(FALSE, 9)) {
  # Start logging
  if (!is.null(log)) sink(log, split=TRUE)
  
  steps <- c("Download tower data", "Download flux data", "Download LiDAR",
             "Calibrate LiDAR constants", "Partition flux",
             "Interpolate meteorology", "Fit AQ curves",
             "Fit Medlyn slopes", "Run energy balance model")
  
  # Initial output
  cat("Output directory:", work_dir, "\n")
  cat("AMF username:", amf_user, "\n")
  cat("AMF email:", amf_email, "\n")
  cat("NEON token is", ifelse(is.null(neon_token), "Missing", "Provided"), "\n")
  
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
  else cat("Skipping tower data download...\n")
  
  #TODO stack tables one at a time, rowbind results, then delete site-months
  # once we are done with them.
  # N.b. this takes about 1.5 hours and requires 60 GB of storage. Warn the 
  # user before proceeding.
  if (!skip[2]) {
    cat("Downloading flux takes ~90 minutes and requires 60 GB of storage. 
        After processing, this is reduced to a 33 MB csv. 
        Are you sure you want to continue?\n")
    cat("1: Continue and download flux\n")
    cat("2: Skip and use \n")
    choice <- readline("Your choice: ")
    if (choice == 1) get_neon_flux(site_meta, neon_token, work_dir)
    else cat("Skipping flux data download...\n")
  } 
  else cat("Skipping flux data download...\n")
  
  if (!skip[3]) get_neon_lidar(site_meta, lidar_dir)
  else cat("Skipping LiDAR data download...\n")
  
  # Fit lidar constants
  if (!skip[4]) fit_neon_lidar_k(site_meta, manual_lai, tower_dir, lidar_dir)
  else cat("Skipping LiDAR constant calibration...\n")
  
  # Partition flux
  raw_flux <- read_csv(file.path(work_dir, "cross_site_flux.csv"))
  if (!skip[5]) partition_neon_flux(site_meta, raw_flux, work_dir)
  else cat("Skipping flux partitioning...\n")
  
  # Within-canopy meteorology
  partition_flux <- read_csv(file.path(work_dir, "cross_site_flux_partition_qc.csv"))
  lad_profiles <- read_csv(file.path(work_dir, "neon_lad_profiles.csv"))
  if (!skip[6]) get_within_canopy_meteorology(site_meta, manual_lai, 
                                              lad_profiles, partition_flux,
                                              work_dir)
  else cat("Skipping meteorology interpolation...\n")
  
  # Fit AQ curves
  if (!skip[7]) fit_aq_curves(site_meta, partition_flux, tower_dir, work_dir)
  else cat("Skipping AQ curve fitting...\n")
  
  
  # Fit Medlyn slopes
  if (!skip[8]) fit_medlyn_slopes(site_meta, manual_lai, partition_flux, 
                                  tower_dir, outdir)
  else cat("Skipping Medlyn slope fitting...\n")
  
  # Run the EB model
  aq_constants <- read_csv(file.path(work_dir, "cross_site_aq_constants.csv"))
  lidar_constants <- read_csv(file.path(work_dir, "neon_lidar_constants.csv"))
  medlyn_constants <- read_csv(file.path(work_dir, "cross_site_medlyn_coefficients.csv"))
  interp_meteo <- read_csv(file.path(work_dir, "cross_site_interpolated_meteorology.csv"))
  
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
  else cat("Skipping energy balance model run...\n")
  
  # Reset logging
  if (!is.null(log)) sink()
}
