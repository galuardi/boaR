% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boa.R
\name{boa}
\alias{boa}
\title{boa}
\usage{
boa(lon, lat, ingrid, nodata = NA, direction = F)
}
\arguments{
\item{lon}{longitude (vector) of the satellite image}

\item{lat}{latitude (vector) of the satellite image}

\item{ingrid}{The satellite data (matrix). If using chlorophyll, transform using log(ingrid))}

\item{nodata}{value representing 'no data'}

\item{direction}{Logical. Should direction be calculated and returned}
}
\value{
either a stand alone asc grid of front gradient data or a list of:
	grdir : ascii grid Gradient direction
   front : Gradient magnitude. In the case of chlorophyll, this is a ratio
}
\description{
boa
}
\details{
Note: These grids are in raster format, as is used in \code{\link{raster}}.
}
\examples{
# none
}
\references{
Belkin, I. M. & O'Reilly, J. E. An algorithm for oceanic front detection in chlorophyll and SST satellite imagery. Journal of Marine Systems, 2009, 78, 319 - 326
}
\keyword{chl,}
\keyword{front,}
\keyword{satellite}
\keyword{sst,}
