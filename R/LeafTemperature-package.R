#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom amerifluxr amf_filter_base
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggrepel geom_label_repel
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @importFrom gridExtra grid.arrange
#' @importFrom gtable gtable_add_grob
#' @importFrom lemon reposition_legend
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom purrr pmap
#' @importFrom sf st_as_sf
#' @importFrom sf st_crs
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom stringr str_split_i
#' @importFrom tigris shift_geometry
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
#' Initially, we used a random sample of NEON DHPs to calculate leaf area index in each of our study sites. Due to large variance in the results, we decided to use the MODIS LAI product instead. These data document our initial attempt with DHPs. `manual_lai_images` documents what the random DHP sample was for reproducibility, and `dhp_lai` records the median LAI we calculated at each site using the `hemispheR` package. See `get_neon_lai.R` for more information.
#' 
NULL

#' @name dhp_lai
#' @rdname manual_lai_images
#' 
NULL

#' @title MODIS leaf area index table
#' @name modis_lai
#' @docType data
#' @description
#' Due to large variation in DHP-derived LAI, we decided to use LAI from MODIS. The process for generating these data was as follows. First, we used AppEARS to derive a point sample of all 8-day 500-m MODIS LAI images from May - September, 2018 - 2022. This request is documented under appears_request in the package directory (or under inst, if you are looking on GitHub). Then, we subset the AppEARS request to those observations where the QC band was zero and calculated site LAI as the 95th percentile of the LAI band. This procedure resulted in better agreement with literature values for LAI than the DHP method, which guided our decision.
NULL
