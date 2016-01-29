\name{boa}
\alias{boa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
boa
}
\description{
Performs the Belkin-O'Reilly smoothing step. This makes use of a sobel filter, land mask and edge correction.
}
\usage{
boa(lon, lat, ingrid, nodata = NA, direction = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{
longitude (vector) of the satellite image
}
  \item{lat}{
latitude (vector) of the satellite image
}
  \item{ingrid}{
The satellite data (matrix). If using chlorophyll, transform using log(ingrid))
}
  \item{nodata}{
value representing 'no data'
}
  \item{direction}{
Logical. Should direction be calculated and returned
}
}
\value{
either a stand alone asc grid of front gradient data or a list of:
\item{grdir}{ascii grid Gradient direction}
\item{front}{Gradient magnitude. In the case of chlorophyll, this is a ratio}

Note: These grids are in asc format, as is used in \code{\link{adehabitat}}. They may be suitable for export to GIS programs but may require some modification in the header (e.g. xll ,yll, cellsize)

}
\references{
Belkin, I. M. & O'Reilly, J. E. An algorithm for oceanic front detection in chlorophyll and SST satellite imagery. Journal of Marine Systems, 2009, 78, 319 - 326
}

\keyword{sst}
\keyword{chl}
\keyword{front}
\keyword{satellite}
