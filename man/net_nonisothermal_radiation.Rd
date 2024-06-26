% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiation.R
\name{net_nonisothermal_radiation}
\alias{net_nonisothermal_radiation}
\title{Net radiation, but allowing for nonisothermal conditions.}
\usage{
net_nonisothermal_radiation(SW_abs, LW_abs_dn, Tl, a_lw = 0.95)
}
\arguments{
\item{SW_abs}{Absorbed downwelling shortwave radiation (W/m2)}

\item{LW_abs_dn}{Absorbed downwelling longwave radiation (W/m2)}

\item{Tl}{Leaf temperature (K)}

\item{a_lw}{Leaf absorptance to longwave radiation (0-1)}

\item{a_sw}{Leaf absorptance to shortwave radiation (0-1)}
}
\value{
Rnet, net radiation (W/m2)
}
\description{
Net radiation, but allowing for nonisothermal conditions.
}
\keyword{internal}
