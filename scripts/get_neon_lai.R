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

calculate_dhp_lai <- function(url, verbose=FALSE, ...) {
  dhp_file1 <- gsub("\\\\", "/", tempfile(fileext=".NEF", tmpdir="data_working"))
  dhp_file2 <- gsub("\\\\", "/", tempfile(fileext=".jpg", tmpdir="data_working"))
  
  curl_command <- make_curl_command(url, dhp_file1)
  dark_command <- make_darktable_command(dhp_file1, dhp_file2)
  
  if (verbose) {
    #cat(curl_command, "\n")
    #cat(dark_command, "\n")
  }
  
  #if (verbose) {cat("Downloading...")}
  system(curl_command, ignore.stdout = TRUE, show.output.on.console = FALSE)
  #if (verbose) {cat("Converting to jpg...\n")}
  system(dark_command, ignore.stdout = TRUE, show.output.on.console = FALSE)
  
  #if (verbose) {cat("Loading image...")}
  dhp <- import_fisheye(dhp_file2,
                        circular=FALSE,
                        display=FALSE,
                        message=FALSE)
  
  #if (verbose) {cat("Binarizing...")}
  dhp_bw <- binarize_fisheye(dhp, method="Otsu",
                             zonal=FALSE, manual=NULL,
                             display=FALSE)
  
  #if (verbose) {cat("Calculate gap fraction...")}
  dhp_gap <- gapfrac_fisheye(dhp_bw, lens="Nikkor-10.5", display=FALSE,
                             nrings=5, nseg=8)
  
  #if (verbose) {cat("Final LAI calculation...\n")}
  dhp_lai <- canopy_fisheye(dhp_gap)
  
  unlink(dhp_file1)
  unlink(dhp_file2)
  
  if (verbose) {
    cat("L =", dhp_lai$L, "\n")
    cat("Le =", dhp_lai$Le, "\n")
  }
  
  return(
    c(dhp_lai, ...)
  )
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
  group_by(siteID, year) %>%
  # Sample at least 10 samples from each group
  slice(sample(min(10, n())))

write_csv(dhp_table, "data_out/neon_lai_metadata_sample.csv")

neon_lai_list <- lapply(
  1:nrow(dhp_table), function(i) {
    cat("[Image", i, "of", nrow(dhp_table), "]\n")
    
    site <- dhp_table$siteID[i]
    year <- dhp_table$year[i]
    url <- dhp_table$imageFileUrl[i]
    
    cat(site, year, "\n")
    cat(url, "\n")
    
    tryCatch(
      calculate_dhp_lai(url, verbose=TRUE, site=site, year=year, srcurl=url),
      error = function(e) list(site=site, year=year, srcurl=url)
    )
  }
)

neon_lai <- bind_rows(neon_lai_list)

write_csv(neon_lai, "data_out/neon_lai.csv")
