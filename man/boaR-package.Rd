\name{boaR-package}
\alias{boaR-package}
\alias{boaR}
\docType{package}
\title{
The Belkin-O'Reilly front detection algorithm.
}
\description{
Belkin-O'Reilly routine for detecting fronts from sea surface temperature and chlorophyll (and probably other) satellite imagery. 
}
\details{
\tabular{ll}{
Package: \tab boaR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-02\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
This package works on satellite imagery imported into R. This package contains three functions which iteratively perfom median filtering and front detection. This is a port from IDL code which operated on flat binary files directly. It is easiest (in my opinion) to import satellite data into R as a netcdf file using the \code{\link{ncdf}} package.   

}
\author{
Benjamin Galuardi

drdrumfish@gmail.com

}
\references{
Belkin, I. M. & O'Reilly, J. E. An algorithm for oceanic front detection in chlorophyll and SST satellite imagery. Journal of Marine Systems, 2009, 78, 319 - 326

}

\examples{
require(fields)

# some good color ramps for gradient mapping
sstgradcol =   colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), interpolate = 'spline', bias = 2.5)
chlgradcol =   colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), interpolate = 'spline', bias = 1.5)

# load some satelite data originally downloaded here (as a netcdf file):  http://coastwatch.pfeg.noaa.gov/erddap/griddap/index.html
data(satdata)

# Blended SST image
attach(sst)
tfront = boa(lon, lat, DATA, direction = T)
par(mfrow=c(1,2))
image.plot(lon, lat, tfront[[1]])
title('Direction')
image.plot(lon, lat, tfront[[2]]) # could use col = sstgradcol(100) for a different scale
title('Gradient Magnitude')
detach(sst)


# SeaWifs chlorophyll image
x11()
attach(swchl)
cfront = boa(swlon, swlat, log(chl), direction = T)  # be sure to take the log of the chl 
par(mfrow=c(1,2))
image.plot(swlon, swlat, cfront[[1]])
title('Direction')
image.plot(swlon, swlat, cfront[[2]], col = chlgradcol(100))
title('Gradient Magnitude (Ratio)')
detach(swchl)

}
