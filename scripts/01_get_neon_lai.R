library(hemispheR)
library(neonUtilities)
library(tidyverse)

fileName <- "data_in/neon_token"
if (file.exists(fileName)) {
  token <- readChar(fileName, file.info(fileName)$size)
} else {
  token <- NA_character_
}

make_curl_command <- function(url, outfile) {
  paste("curl", url, "-L --output", outfile)
}

make_darktable_command <- function(input, output) {
  paste("bash -c ", "'darktable-cli ", input, " ", output, "'", sep='')
}

download_convert_dhp <- function(url, verbose=FALSE, ...) {
  dhp_file1 <- gsub("\\\\", "/", tempfile(fileext=".NEF", tmpdir="data_working"))
  dhp_file2 <- gsub("\\\\", "/", tempfile(fileext=".jpg", tmpdir="data_working"))
  
  curl_command <- make_curl_command(url, dhp_file1)
  dark_command <- make_darktable_command(dhp_file1, dhp_file2)
  
  if (verbose) {cat("Downloading...\n", curl_command, "\n")}
  system(curl_command, ignore.stdout = TRUE, show.output.on.console = FALSE)
  if (verbose) {cat("Converting to jpg...\n", dark_command, "\n")}
  system(dark_command, ignore.stdout = FALSE, show.output.on.console = TRUE)
  
  unlink(dhp_file1)
  
  return(dhp_file2)
}

calculate_binarized_lai <- function(dhp_bw, lens="Nikkor-10.5",
                                    display=FALSE, nrings=5, nseg=8, ...) {
  
  dhp_gap <- gapfrac_fisheye(dhp_bw, lens="Nikkor-10.5", display=FALSE,
                             nrings=5, nseg=8)
  
  dhp_lai <- canopy_fisheye(dhp_gap)
  
  return(c(dhp_lai, ...))
}

calculate_dhp_lai_auto <- function(url, verbose=FALSE, ...) {
  
  dhp_file <- download_convert_dhp(url, verbose=verbose)
  
  #if (verbose) {cat("Loading image...")}
  dhp <- import_fisheye(dhp_file,
                        circular=FALSE,
                        display=FALSE,
                        message=FALSE)
  
  #if (verbose) {cat("Binarizing...")}
  dhp_bw <- binarize_fisheye(dhp, method="Otsu",
                             zonal=FALSE, manual=NULL,
                             display=FALSE)
  
  
  result <- calculate_binarized_lai(dhp_bw, ...)
  
  unlink(dhp_file)
  
  return(result)
} 

calculate_dhp_lai_manual <- function(url, verbose=FALSE, ...) {
  dhp_file <- download_convert_dhp(url, verbose=verbose)
  
  dhp <- import_fisheye(dhp_file,
                        circular=FALSE,
                        display=FALSE,
                        message=FALSE)
  
  thresh <- as.numeric(autothresholdr::auto_thresh(
    round(as.matrix(dhp)), method="Otsu"
  ))
  dhp_bw <- binarize_fisheye(dhp, manual=thresh, display=FALSE)
  par(mfrow=c(1, 2))
  while (TRUE) {
    cat("Current thresh:", thresh, "\n")
    dhp_bw <- binarize_fisheye(dhp, manual=thresh, display=FALSE)
    
    plot(dhp, col=grey.colors(50), main="Original")
    plot(dhp_bw, col=grey.colors(2), main="Thresholded")
    
    inp <- readline("[A]ccept, [R]eject, or input new value: ")
    
    if (inp == "A") {break}
    if (inp == "R") {return(list(...))}
    
    newinp <- suppressWarnings(as.numeric(inp))
    
    if (!is.na(newinp)) {
      thresh <- newinp
    }
  }
  par(mfrow=c(1, 1))
  
  result <- calculate_binarized_lai(dhp_bw, ...)
  
  unlink(dhp_file)
  
  return(result)
}

sites <- read_csv("data_working/neon_site_metadata.csv") %>% pull(site_neon)

dhp_product <- loadByProduct(
  "DP1.10017.001",
  site=sites,
  token=token,
  tabl="dhp_perimagefile",
  package="basic",
  check.size=FALSE
)

dhp_table <- dhp_product$dhp_perimagefile %>%
  filter(imageType == "overstory",
         cameraOrientation == "upward") %>%
  mutate(year = year(startDate)) %>%
  filter(year >= 2019) %>%
  group_by(siteID) %>%
  # Sample at least 10 images from each site
  slice(sample(min(10, n())))

write_csv(dhp_table, "data_out/neon_sampled_dhp_images.csv")

neon_lai_list <- lapply(
  1:nrow(dhp_table), function(i) {
    cat("[Image", i, "of", nrow(dhp_table), "]\n")
    
    site <- dhp_table$siteID[i]
    year <- dhp_table$year[i]
    url <- dhp_table$imageFileUrl[i]
    
    cat(site, year, "\n")
    cat(url, "\n")
    
    tryCatch(
      calculate_dhp_lai_manual(url, verbose=TRUE, site=site, year=year, srcurl=url),
      error = function(e) list(site=site, year=year, srcurl=url)
    )
  }
)

neon_lai_list_filter <- neon_lai_list[lapply(neon_lai_list, class) == "list"]
neon_lai <- bind_rows(neon_lai_list_filter)

median_lai <- neon_lai %>% group_by(site) %>% summarize(L_dhp = median(L))

# Compare the DHP-derived LAI with other estimates in the literature
known_lais <- read_csv("data_working/neon_site_metadata.csv") %>%
  select(site_neon, lai) %>% filter(!is.na(lai)) %>%
  left_join(median_lai, by=c("site_neon"="site"))

correction_factor <- known_lais %>%
  # Geometric mean error ratio
  mutate(error_ratio = lai / L_dhp) %>%
  pull(error_ratio) %>%
  log() %>% mean() %>% exp()

neon_lai_corrected <- median_lai %>%
  # Very sparse ecosystems give an LAI < 1, but a value this low doesn't really
  # make sense for the energy balance model. Consider the LAI = 1 system to be
  # a section of forest where there are enough trees to yield a meaningful LAI.
  mutate(L_dhp_corrected = pmax(1, L_dhp * correction_factor))

write_csv(neon_lai_corrected, "data_out/neon_sampled_dhp_lai.csv")
