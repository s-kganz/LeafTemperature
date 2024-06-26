% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leaf.R
\name{latent_heat_imposed}
\alias{latent_heat_imposed}
\title{Imposed latent heat flux. This is the latent heat flux when air temperature
and leaf temperature are well-coupled.}
\usage{
latent_heat_imposed(gtot, vpd, Pa, l = eb_constants_$l)
}
\arguments{
\item{gtot}{Total conductance (mol / m2 sec)}

\item{vpd}{Vapor pressure deficit (kPa)}

\item{Pa}{Air pressure (kPa)}

\item{l}{Heat of vaporization of water (J / mol)}
}
\value{
LE_imp, stomatally-imposed latent heat flux.
}
\description{
Imposed latent heat flux. This is the latent heat flux when air temperature
and leaf temperature are well-coupled.
}
\keyword{internal}
