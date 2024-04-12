#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom amerifluxr amf_filter_base
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @importFrom lemon reposition_legend
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom purrr pmap
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom stringr str_split_i
## usethis namespace: end
NULL

#' @title Site 
#' @name site_meta
#' @docType data
#' @description
#' This is a table of site attributes used across many of the functions in this package.
#' 
NULL

#' @title Site leaf area index tables
#' @name manual_lai_images
#' @docType data
#' @description
#' We used a random sample of NEON DHPs to calculate leaf area index in each of our study sites. `manual_lai_images` documents what that random sample was for reproducibility, and `manual_lai` records the median LAI we calculated at each site using the `hemispheR` package. See `get_neon_lai.R` for more information.
#' 
NULL

#' @name manual_lai
#' @rdname manual_lai_images
#' 
NULL
