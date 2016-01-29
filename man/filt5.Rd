\name{filt5}
\alias{filt5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filt5
}
\description{
Finds maximum of a 5x5 sliding window. If the central pixel is the maximum, this is flagged as a one. All other pixels are flagged as zero. Usually this function is not used directly.
}
\usage{
filt5(lon, lat, ingrid, nodata = NA)
}

\arguments{
  \item{lon}{
longitude (vector) of the satellite image
}
  \item{lat}{
latitude (vector) of the satellite image
}
  \item{ingrid}{
The satellite data (matrix)
}
  \item{nodata}{
value representing 'no data'
}
}
\value{
returns a grid of zeros and ones
}
\seealso{
\code{\link{boa}}
}
