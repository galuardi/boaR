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

