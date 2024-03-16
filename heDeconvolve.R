heDeconvolve <- function(
    
  ## Input
  imArray,
  
  HStain =  c(0.65, 0.704, 0.286), #'#594BB6' upon 1-Hx
  EStain = c(0.072, 0.990, 0.105) #'#ED03E4' upon 1-Ey
){
  
  if(length(dim(imArray)) != 3){
    cat('\n\nERROR: EXPECTING AN ARRAY WITH 3 DIMENSIONS!')
    break
  }
  if(dim(imArray)[3] != 3){
    cat('\n\nERROR: EXPECTING 3 COLOURS IN 3 DIMENSION!')
    break
  }
  if(length(HStain) != 3 | length(EStain) !=3 ){
    cat("\n\nERROR: EXPECTING H AND/OR E STAINS TO BE A RGB VECTOR!")
    break
  }
  require(MASS)
  v <- matrix(c(HStain, EStain), byrow = FALSE, ncol=2, nrow=3)
  im <- imArray
  od <- do.call(cbind, lapply(1:3, function(i) as.vector(1-im[,,i]) - mean(1-im[,,i]) )) #Subtract mean intensity to ensure both coefficients centred at 0
  s <- MASS::ginv(t(v) %*% v) %*% t(v) %*% t(od)
  imh <- matrix(s[1,], nrow=nrow(im), ncol=ncol(im))
  ime <- matrix(s[2,], nrow=nrow(im), ncol=ncol(im))
  imh[imh<0] <- 0
  ime[ime<0] <- 0
  imh = imh / max(imh)
  ime = ime / max(ime)
  
  h <- -Reduce('+', lapply(heList, function(x){ 
    
    x[x<0] <- 0
    x <- x / max(x)
    y <- x * log2(x)
    y[x==0] <- 0
    return(y)
    
  }))
  
  OUT <- list(
    'H' = imh,
    'E' = ime,
    'ENTROPY' = h
  )
  
  return(OUT)
}