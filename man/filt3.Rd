\name{filt3}
\alias{filt3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filt3
}
\description{
Performs a median filter on a grid of satelite data using a 3x3 sliding window. This is used following \code{\link{filt5}} where pixels of maximum value were found over a 5x5 sliding window. these maximum pixel values are  retained in this step. Ususally, this is not called directly.
}
\usage{
filt3(lon, lat, ingrid, grid5 = grid5(lon, lat, ingrid))
}

\arguments{
  \item{lon}{
longitude (vector) of the satellite image
}
  \item{lat}{
latitude (vector) of the satellite image
}
  \item{grid5}{
results of \code{\link{filt5}}
}
}
\value{
returns a median smoothed grid of satellite data
}
\seealso{
\code{\link{boa}}
}