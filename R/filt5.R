#' filt5
#' @description  Finds maximum of a 5x5 sliding window. If the central pixel is the maximum, this is flagged as a one. All other pixels are flagged as zero. Usually this function is not used directly.
#' @usage filt5(lon, lat, ingrid, nodata = NA)
#' @param lon longitude (vector) of the satellite image
#' @param lat latitude (vector) of the satellite image
#' @param ingrid The satellite data (matrix)
#' @param nodata value representing 'no data'
#'
#' @return returns a grid of zeros and ones
#' @seealso \code{\link{boa}}
#' @export
#'
#' @examples 
#' # none
filt5 <-
function(lon , lat, ingrid, nodata = NA){
#======================================================#
# Find peaks in 5 directions. FLag as 5
#======================================================#
 #outgrid = matrix(0,nrow=dim(ingrid)[1],ncol=dim(ingrid)[2] ) # make all zeros..
 nodatidx = ingrid[ingrid==nodata]
 ingrid[nodatidx] = -9999
 outgrid = ingrid*0 # make all zeros..
 l1 = length(lon)
 l2 = length(lat)
for(i in 3:(l1-2)){
 for(j in 3:(l2-2)){
   subg = ingrid[(i-2):(i+2),(j-2):(j+2)]
if(sum(is.na(subg))==25){
outgrid[i,j] = 0
}else{
vec = as.vector(subg)
ma = which.max(subg)
mi = which.min(subg)

   if(ma==13||mi==13){
 outgrid[i,j] = 1
}else{
  outgrid[i,j] = 0
}
}
   }
  }
outgrid[nodatidx] = nodata  
outgrid
}

