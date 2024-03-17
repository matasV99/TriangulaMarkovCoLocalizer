fftshift <- function(img_ff, dim = -1) {
  
  rows <- dim(img_ff)[1]    
  cols <- dim(img_ff)[2]    
  
  swap_up_down <- function(img_ff) {
    rows_half <- ceiling(rows/2)
    return(rbind(img_ff[((rows_half+1):rows), (1:cols)], img_ff[(1:rows_half), (1:cols)]))
  }
  
  swap_left_right <- function(img_ff) {
    cols_half <- ceiling(cols/2)
    return(cbind(img_ff[1:rows, ((cols_half+1):cols)], img_ff[1:rows, 1:cols_half]))
  }
  
  if (dim == -1) {
    img_ff <- swap_up_down(img_ff)
    return(swap_left_right(img_ff))
  }
  else if (dim == 1) {
    return(swap_up_down(img_ff))
  }
  else if (dim == 2) {
    return(swap_left_right(img_ff))
  }
  else {
    stop("Invalid dimension parameter")
  }
}

# cat('Inverse shift pixels for FFT analysis...')
ifftshift <- function(img_ff, dim = -1) {
  
  rows <- dim(img_ff)[1]    
  cols <- dim(img_ff)[2]    
  
  swap_up_down <- function(img_ff) {
    rows_half <- floor(rows/2)
    return(rbind(img_ff[((rows_half+1):rows), (1:cols)], img_ff[(1:rows_half), (1:cols)]))
  }
  
  swap_left_right <- function(img_ff) {
    cols_half <- floor(cols/2)
    return(cbind(img_ff[1:rows, ((cols_half+1):cols)], img_ff[1:rows, 1:cols_half]))
  }
  
  if (dim == -1) {
    img_ff <- swap_left_right(img_ff)
    return(swap_up_down(img_ff))
  }
  else if (dim == 1) {
    return(swap_up_down(img_ff))
  }
  else if (dim == 2) {
    return(swap_left_right(img_ff))
  }
  else {
    stop("Invalid dimension parameter")
  }
}

tileFFTImage <- function(
    
  ## INPUT
  imArray,
  
  tileSize = 60,
  overlapSize = 30,
  verbose = T
  
){
  
  ### SET UP
  if(is.list(imArray)){
    imList <- imArray
    dim_check = T
    for(i in 1:length(imList)){
      dim_check = dim_check & all(
        dim(imList[[1]]) == dim(imList[[i]])
      )
    }
    if(!dim_check){
      cat('\n\nERROR: IMAGE LIST PROVIDED DOES NOT CONTAIN MATRICES OF EQUAL DIMENSIONS!')
      break
    }
  }else if(length(dim(imArray))==2){
    imList <- list(imArray)
  }else if(length(dim(imArray))==3){
    imList <- lapply(1:dim(imArray)[3], function(i){ return(imArray[,,i]) })
  }else{
    cat('\n\nERROR: INVALID IMAGE INPUT!')
    break
  }
  dim_m = dim(imList[[1]])
  if(tileSize < 1 | overlapSize <0){
    cat('\n\nERROR: INVALID tileSize AND/OR overlapSize!')
    break
  }
  step_size = tileSize - overlapSize
  if(step_size <0){
    cat('\n\nERROR: overlapSize CANNOT BE BIGGER THAN tileSize!')
    break
  }
  require(abind)
  require(RcppML)
  ###
  
  ### TILING
  if(verbose){cat('\nGetting tile coordinates...')}
  cutX <- as.integer(cut(1:dim_m[1], seq(0, dim_m[1], by=step_size)+1, include.lowest = T))
  cutX[1] <- 0
  cutX[is.na(cutX)] <- max(cutX[!is.na(cutX)]) + 1
  cutY <- as.integer(cut(1:dim_m[2], seq(0, dim_m[2], by=step_size)+1, include.lowest = T))
  cutY[1] <- 0
  cutY[is.na(cutY)] <- max(cutY[!is.na(cutY)]) + 1
  
  if(verbose){cat('\nGetting Hanning window...')}
  M = tileSize
  n = seq(1-M, M, by=2)
  col_hanning <- row_hanning <- 0.5 + 0.5 * cos(pi * n/(M-1))
  hanning <- matrix(row_hanning, ncol = 1) %*% matrix(col_hanning, nrow = 1)
  
  OUT <- NULL
  for(d in 1:length(imList)){
    if(verbose){ cat(paste0('\n', d, ' of ', length(imList), '...\n'))}
    imRast <- imList[[d]]
    idx = 1
    RES <- array(0, dim=c(max(which(tabulate(cutX)==step_size)), max(which(tabulate(cutY)==step_size)), tileSize^2))
    for(x in unique(cutX)){
      if(sum(cutX==x) != step_size){next}
      if(verbose){ cat(paste0(x, ' of ', max(cutX), '...'))}
      xstart = min(which(cutX==x))
      xend = xstart + tileSize-1
      if(xend > length(cutX)){next}
      for(y in unique(cutY)){
        if(sum(cutY==y) != step_size){next}
        ystart = min(which(cutY==y))
        yend = ystart + tileSize -1
        if(yend > length(cutY)){next}
        tile <- imRast[xstart:xend, ystart:yend]
        tile_fft <- fftshift(fft(tile * hanning))
        RES[x,y,] <- Mod(tile_fft)
      }
    }
    if(is.null(OUT)){
      OUT <- RES
    }else{
      OUT <- abind(OUT, RES, along = 3)
    }
  }
  ###
  
  return(OUT)
}