boa <-
function(lon, lat, ingrid, nodata = NA, direction = F){
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
land = naidx*1
land[land==1] = NaN
land[!is.nan(land)] = 1
# use adehabitat library for image erosion (expand land pixels slightly)
mask = morphology(as.asc(land, min(lon), yll = min(lat), cellsize = diff(lat)[1]), 'erode', 1)
# remove edge pixels
l1 = length(lon)
l2 = length(lat)
midx = mask*NaN
midx[5:(l1-3), 5:(l2-3)]=1
mask = mask*midx
# account for the edge effect
front = front *mask

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
list(grdir = GRAD_DIR*mask, front = front)
}
else{
front
}
}

