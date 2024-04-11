starting <- function(msg) {
  message("Starting ", msg)
}

skipping <- function(msg) {
  message("Skipping ", msg)
}

#' Re-run leaf temperature modeling pipeline
#'
#' @param amf_user Ameriflux username for data download
#' @param amf_email Ameriflux email for data download
#' @param neon_token Token for NEON API. Can be NA, but this will makes your downloads take longer.
#' @param work_dir Working directory to write tables.
#' @param log If not `NULL`, all output will be written to this file.
#' @param skip 11-character string describing whether you want to skip certain analysis steps. This can be useful if you hit an error halfway through and want to skip all the steps that succeeded initially.
#'
#' @details
#' This function re-runs the leaf temperature modeling pipeline we describe in the manuscript. It begins by downloading ~2 GB of data to your machine from Ameriflux and NEON. If you run the whole thing, you should expect it to take 2-3 hours. Most of that time is just downloading data, the actual modeling part takes just a few minutes. We have done our best to make the code robust, but ti is always possible that the analysis breaks halfway through. If that happens, you can adjust the `skip` parameter to ignore the steps that have already completed.
#'
#' The skip parameter is an 11-character string. If the i-th character is 'y', then the i-th analysis step is skipped. If a step is skipped, it is assumed that the relevant data has already been written to `work_dir` with the expected filename. The order of steps is as follows:
#'  \enumerate{
#'   \item [get_amf_tower_data()]
#'   \item [get_neon_flux()]
#'   \item [get_neon_lidar()]
#'   \item [fit_neon_lidar_k()]
#'   \item [partition_neon_flux()]
#'   \item [get_within_canopy_meteorology()]
#'   \item [fit_aq_curves()]
#'   \item [fit_medlyn_slopes()]
#'   \item [run_neon_energy_balance()]
#'   \item [model_radiometer_comparison()]
#'   \item [gs_gbh_sensitivity()]
#'  }
#'
#'  See the linked functions for more information. Once the analysis finishes, you can run [write_all_figures()] to generate the figures we show in the manuscript.
#'
#' @export
#'
run_analysis <-
  function(amf_user,
           amf_email,
           neon_token,
           work_dir,
           log = NULL,
           skip = "nnnnnnnnnnn") {
    # Start logging
    if (!is.null(log))
      sink(log, split = TRUE)
    
    steps <-
      c(
        "Download tower data",
        "Download flux data",
        "Download LiDAR",
        "Calibrate LiDAR constants",
        "Partition flux",
        "Interpolate meteorology",
        "Fit AQ curves",
        "Fit Medlyn slopes",
        "Run energy balance model",
        "Model vs. radiometer comparison",
        "Conductance sensivity analysis"
      )
    
    skip <- (strsplit(skip, "")[[1]]) == "y"
    
    # Initial output
    message("Output directory: ", work_dir)
    message("AMF username: ", amf_user)
    message("AMF email: ", amf_email)
    message("NEON token is ",
            ifelse(
              is.null(neon_token) | is.na(neon_token),
              "Missing",
              "Provided"
            ))
    
    # Set up directories
    if (!dir.exists(work_dir))
      dir.create(work_dir)
    tower_dir <- file.path(work_dir, "amf_tower_data")
    lidar_dir <- file.path(work_dir, "neon_neartower_lidar")
    model_dir <- file.path(work_dir, "model_run")
    
    dir.create(tower_dir)
    dir.create(lidar_dir)
    dir.create(model_dir)
    
    # Load site metadata and LAI. Technically this isn't necessary but just FYI
    # these are included with the package.
    data("site_meta")
    data("manual_lai")
    
    # Download all the data we need
    if (!skip[1]) {
      starting(steps[1])
      get_amf_tower_data(site_meta, work_dir, amf_user, amf_email)
    }
    else
      skipping(steps[1])
    
    if (!skip[2]) {
      starting(steps[2])
      get_neon_flux(site_meta, neon_token, work_dir)
    }
    else
      skipping(steps[2])
    
    if (!skip[3]) {
      starting(steps[3])
      get_neon_lidar(site_meta, lidar_dir)
    }
    else
      skipping(steps[3])
    
    # Fit lidar constants
    if (!skip[4]) {
      starting(steps[4])
      fit_neon_lidar_k(site_meta, manual_lai, tower_dir, lidar_dir, work_dir)
    }
    else
      skipping(steps[4])
    
    # Partition flux
    raw_flux <- read_csv(file.path(work_dir, "cross_site_flux.csv"),
                         show_col_types = FALSE)
    if (!skip[5]) {
      starting(steps[5])
      partition_neon_flux(site_meta, raw_flux, tower_dir, work_dir)
    }
    else
      skipping(steps[5])
    
    # Within-canopy meteorology
    partition_flux <-
      read_csv(file.path(work_dir, "cross_site_flux_partition_qc.csv"),
               show_col_types = FALSE)
    lad_profiles <-
      read_csv(file.path(work_dir, "neon_lad_profiles.csv"),
               show_col_types = FALSE)
    if (!skip[6]) {
      starting(steps[6])
      get_within_canopy_meteorology(site_meta,
                                    manual_lai,
                                    lad_profiles,
                                    partition_flux,
                                    work_dir)
    }
      else
        skipping(steps[6])
      
      # Fit AQ curves
      if (!skip[7]) {
        starting(steps[7])
        fit_aq_curves(site_meta, partition_flux, tower_dir, work_dir)
      }
      else
        skipping(steps[7])
      
      
      # Fit Medlyn slopes
      if (!skip[8]) {
        starting(steps[8])
        fit_medlyn_slopes(site_meta, manual_lai, partition_flux,
                          tower_dir, work_dir)
      }
      else
        skipping(steps[8])
      
      # Run the EB model
      aq_constants <-
        read_csv(file.path(work_dir, "cross_site_aq_constants.csv"),
                 show_col_types = FALSE)
      lidar_constants <-
        read_csv(file.path(work_dir, "neon_lidar_constants.csv"),
                 show_col_types = FALSE)
      medlyn_constants <-
        read_csv(file.path(work_dir, "cross_site_medlyn_coefficients.csv"),
                 show_col_types = FALSE)
      interp_meteo <-
        read_csv(
          file.path(work_dir, "cross_site_interpolated_meteorology.csv"),
          show_col_types = FALSE
        )
      
      if (!skip[9]) {
        starting(steps[9])
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
      }
      else
        skipping(steps[9])
      
      eb_result <- read_csv(file.path(model_dir, "cross_site_eb.csv"),
                            show_col_types = FALSE)
      rad_tcan <- read_csv(file.path(work_dir, "cross_site_tcan.csv"),
                           show_col_types = FALSE)
      
      if (!skip[10]) {
        starting(steps[10])
        model_radiometer_comparison(eb_result, rad_tcan, model_dir)
      }
      else
        skipping(steps[10])
      
      if (!skip[11]) {
        starting(steps[11])
        gs_gbh_sensitivity(model_dir)
      }
      else
        skipping(steps[11])
      
      # Reset logging
      if (!is.null(log)) {
        sink()
      }
    }
    