#' filt3
#'
#' @param lon longitude (vector) of the satellite image
#' @param lat latitude (vector) of the satellite image
#' @param ingrid The satellite data (matrix)
#' @param grid5 results of \code{\link{filt5}}
#'
#' @return returns a median smoothed grid of satellite data
#' @seealso \code{\link{boa}}
#' @usage filt3(lon, lat, ingrid, grid5 = grid5(lon, lat, ingrid))
#' @export
#'
#' @examples
#' # none here
filt3 <-
function(lon , lat, ingrid, grid5){
#======================================================#
# Find peaks in 3 directions. FLag as 3
#======================================================#
 outgrid = ingrid*0 # make all zeros..
 l1 = length(lon)
 l2 = length(lat)
for(i in 3:(l1-2)){
 for(j in 3:(l2-2)){
  if((grid5[i,j]==0)==T){
 subg = ingrid[(i-1):(i+1),(j-1):(j+1)]
if(sum(is.na(subg))==9){
outgrid[i,j] = ingrid[i,j]
}else{
vec = as.vector(subg)
ma = which.max(subg)
mi = which.min(subg)

   if((ma==5||mi==5)==T){
 outgrid[i,j] = median(subg,na.rm=T)  # apply median filter if there is peak 3
}else{
  outgrid[i,j] = ingrid[i,j]
}
}
}else{
 outgrid[i,j] = ingrid[i,j]
}

 }
}
outgrid
}

