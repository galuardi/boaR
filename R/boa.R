#' boa
#'
#' @param lon longitude (vector) of the satellite image
#' @param lat latitude (vector) of the satellite image
#' @param ingrid The satellite data (matrix). If using chlorophyll, transform using log(ingrid))
#' @param nodata value representing 'no data'
#' @param direction Logical. Should direction be calculated and returned
#'
#' @return either a stand alone asc grid of front gradient data or a list of:
#'		grdir : ascii grid Gradient direction
#'    front : Gradient magnitude. In the case of chlorophyll, this is a ratio
#' @references Belkin, I. M. & O'Reilly, J. E. An algorithm for oceanic front detection in chlorophyll and SST satellite imagery. Journal of Marine Systems, 2009, 78, 319 - 326
#' @details Note: These grids are in raster format, as is used in \code{\link{raster}}. 
#' @keywords sst, chl, front, satellite
#' @export
#'
#' @examples # none
boa <-
function(lon, lat, ingrid, nodata = NA, direction = F){

## Workhorse filter from EBImage. Modified so we don't need colorspace and other annoying requirements.		
	.filter2 <-
		function (x, filter) 
		{
			# Workhorse filter from EBImage. Modified so we don't need colorspace and other annoying requirements.
			validObject(x)
			validObject(filter)
			#if (colorMode(x) == TrueColor) 
			#   stop("this method doesn't support the 'TrueColor' color mode")
			dx = dim(x)
			#cmx = colorMode(x)
			df = dim(filter)
			if (any(df%%2 == 0)) 
				stop("dimensions of 'filter' matrix must be odd")
			if (any(dx[1:2] < df)) 
				stop("dimensions of 'x' must be bigger than 'filter'")
			cx = dx%/%2
			cf = df%/%2
			wf = matrix(0, nr = dx[1], nc = dx[2])
			wf[(cx[1] - cf[1]):(cx[1] + cf[1]), (cx[2] - cf[2]):(cx[2] + 
																													 	cf[2])] = filter
			wf = fft(wf)
			dim(x) = c(dx[1:2], prod(dx)/prod(dx[1:2]))
			index1 = c(cx[1]:dx[1], 1:(cx[1] - 1))
			index2 = c(cx[2]:dx[2], 1:(cx[2] - 1))
			pdx = prod(dim(x)[1:2])
			y = apply(x, 3, function(xx) {
				dim(xx) = dx[1:2]
				Re(fft(fft(xx) * wf, inverse = T)/pdx)[index1, index2]
			})
			dim(y) = dx
			y
		}	
	
	# browser()
#======================================================#
# Main BOA algorithm
#======================================================#
# require(adehabitat)
# x component (longitudinal kernel)
gx = matrix(c(-1,0,1,-2,0,2,-1,0,1),nrow=3, byrow=T)
# y component (latitudinal kernel)
gy = matrix(c(1,2,1,0,0,0,-1,-2,-1),nrow=3, byrow=T)
# filt 5 and 35 don't like NA's... but final steps are ok with it!

# if nodata is numeric, this will take care of it..
ingrid[ingrid == nodata] = -9999
# if something else: 
if(any(is.infinite(ingrid)|is.nan(ingrid)|is.na(ingrid))){
ingrid[is.na(ingrid)] = -9999
ingrid[is.nan(ingrid)] = -9999
ingrid[is.infinite(ingrid)] = -9999
}

# do the median filtering
grid5 = filt5(lon, lat, ingrid, nodata = nodata)
grid35 = filt3(lon, lat, ingrid, grid5)

# make an index of bad values and land pixels.
grid35[grid35==-9999] = NA
naidx = is.na(grid35)
# convert these to zeros for smoothing purposes
grid35[naidx]=0
# perform the smoothing (Sobel filter)
tgx = .filter2(grid35, gx)
tgy = .filter2(grid35, gy)

# NOTE!! IDL CONVOL uses NORMALIZE function which defaults to the abs(sum(KERNEL)), in this case gx and gy
tx = tgx/sum(abs(as.vector(gx)), na.rm=T)
ty = tgy/sum(abs(as.vector(gy)), na.rm=T)
front = (sqrt(tx^2+ty^2))

#======================================================#
# landmask and edge dilation
#======================================================#
land = naidx*1
# land = naidx*1
land[land==1] = NaN
land[!is.nan(land)] = 1
# # use adehabitat library for image erosion (expand land pixels slightly)
# mask = morphology(as.asc(land, min(lon), yll = min(lat), cellsize = diff(lat)[1]), 'erode', 1)
# # remove edge pixels
# l1 = length(lon)
# l2 = length(lat)
# midx = mask*NaN
# midx[5:(l1-3), 5:(l2-3)]=1
# mask = mask*midx
# # account for the edge effect
# front = front *mask

#======================================================#
# landmask and edge dilation using raster!
#======================================================#


l1 = length(lon)
l2 = length(lat)

midx = land*NaN

midx[5:(l1-3), 5:(l2-3)] = 1

# land = land*midx

land = (land*midx)

# ssr = flip(t(raster(front)), 2)
ssr = flip(raster((t(front))),2)
extent(ssr) = c(min(lon), max(lon), min(lat), max(lat))

mask = flip(focal(raster(t(land)), w = matrix(c(0,0,0,0,1,0,0,0,0), nrow=3, ncol=3)), direction = 2)
extent(mask) = c(min(lon), max(lon), min(lat), max(lat))

front = mask*ssr


if(direction==T){
# ;   ************************************
# ;   *** Calculate Gradient Direction ***
# ;   ************************************
    GRAD_DIR = atan2(tgy@.Data, tgx@.Data)

# ;===> change radians to degrees
GRAD_DIR = GRAD_DIR*180/pi

# ;===> Adjust to 0-360 scheme (make negative degrees positive)
OK = which(GRAD_DIR < 0)
if(length(OK)>1)GRAD_DIR[OK] = 360 - abs(GRAD_DIR[OK])

# ;===> Convert degrees so that 0 degrees is North and East is 90 degrees
GRAD_DIR = (360 - GRAD_DIR + 90) %% 360

GRAD_DIR = flip(raster(t(GRAD_DIR)), direction = 2)
extent(GRAD_DIR) = extent(front)

list(grdir = GRAD_DIR*mask, front = front)
}
else{
front
}
}

