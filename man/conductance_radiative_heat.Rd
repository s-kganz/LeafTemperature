% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leaf.R
\name{conductance_radiative_heat}
\alias{conductance_radiative_heat}
\title{Leaf conductance to radiative heat transfer. Note that the output is in
molar units here for compatability with other functions.}
\usage{
conductance_radiative_heat(
  Ta,
  rho,
  rho_mol,
  a_lw = 0.95,
  cp = eb_constants_$cp,
  sb = eb_constants_$sb
)
}
\arguments{
\item{Ta}{Air temperature (K)}

\item{rho}{Mass density of dry air (kg / m^3)}

\item{rho_mol}{Molar density of dry air (mol / m^3)}

\item{a_lw}{Longwave absorptance (-)}

\item{cp}{Heat capacity of dry air (J / kg K)}

\item{sb}{Stefan-Boltzmann constant (W m-2 K-4)}
}
\value{
gR, leaf conductance to radiative heat transfer (mol / m2 sec)
}
\description{
Leaf conductance to radiative heat transfer. Note that the output is in
molar units here for compatability with other functions.
}
\keyword{internal}
