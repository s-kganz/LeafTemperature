
# Modeling leaf temperature in conifer forests.

This package contains analysis scripts and data described in our manuscript on modeling leaf energy balance in conifer forests. The goal of this package is to make it so that (with sufficient patience) you can reproduce the figures in the manuscript with minimal coding. We hope that this will encourage others to verify our results, extend the model to other study sites, or just learn a little more about leaf energy budgets.

## Installation

For various reasons, this package will not be submitted to CRAN. Instead, please install it directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("s-kganz/LeafTemperature")
```

## Example

Once you have the package installed, reproducing our analysis is straightforward.

``` r
library(LeafTemperature)

run_analysis(
  # Ameriflux API credentials
  your_amf_username,
  your_amf_email,
  # NEON API credentials
  your_neon_token,
  # Scratch directory to save tables
  "scratch"
)

write_all_figures(
  site_meta,
  # Working directory used in run_analysis()
  "scratch",
  # Output location for figures
  "scratch/figures"
)
```

And that's it! Be aware that the analysis pipeline downloads raw data from Ameriflux and NEON, partitions eddy covariance flux, calibrates various constants, and then runs the energy balance model itself. This can take ~ 1 hour to run and, upon completion, the scratch directory will have ~ 2 GB in it.

Also, due to platform-dependent issues with deleting files, flux download just pulls a processed CSV from a Google Cloud bucket. The original function we used to make the table is provided, but it was too temperamental to get working on different machines.

