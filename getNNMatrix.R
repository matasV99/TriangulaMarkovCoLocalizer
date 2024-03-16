getNNMatrix <- function( 
    
  ## INPUT
  x, ## X coordinates
  y, ## Y coordinates
  cell_label, ## Cell type labels
  
  ## Delaunay Triangulation Option: fast, and number of times run depends on delaunayNNDegrees
  delaunayTriangulation = T,
  delaunayTriangulationDistanceThreshold = 20,
  delaunayNNDegrees = c(1),
  
  ## Distances option: much slower, but only done once
  euclideanDistances = NULL,
  
  ## For tripack::tri.mesh function
  duplicate = "error", 
  jitter = 10^-12, 
  jitter.iter = 6, 
  jitter.random = FALSE,
  
  verbose = TRUE
){
  
  
  ###
  if(length(x) != length(y) | length(x) != length(cell_label) ){
    cat('\n\nERROR: X, Y, AND CELL_LABELS REQUIRE THE SAME NUMBER OF ENTRIES!')
    break
  }
  
  if(!delaunayTriangulation & is.null(euclideanDistances)){
    cat("\n\nERROR: EITHER SPECIFY DISTANCES TO SEARCH FOR NEAREST NEIGHBOURS, OR SET delaunayTriangulation TO TRUE!")
    break
  }
  ###
  
  
  ###
  CLB <- as.integer(factor(cell_label))
  OUT <- list()
  all_coords <- x + (y * 1i)
  # idx <- x + max(x) * (y-min(y))
  ###
  
  ###
  if(delaunayTriangulation){
    
    if(verbose){ cat('\n\nRunning delaunay triangulation approach for getting nearest neighbours...') }
    require(tripack)
    dt <- tripack::tri.mesh(
      x=x, y=y, duplicate = duplicate, jitter = jitter, jitter.iter = jitter.iter, jitter.random = jitter.random
    )
    nns_raw <- tripack::neighbours(dt)
    convex_hull <- tripack::convex.hull(dt)
    
    if(verbose){ cat('\nTabulating...') }
    
    DTS <- list()
    for(j in 1:length(delaunayNNDegrees)){
      DEGREE <- delaunayNNDegrees[j]
      cat(paste0('\nExamining neighbours of ', DEGREE, ' degree...'))
      
      ## Iteratively add more neighbours
      di = 1
      while(di <= DEGREE){
        if(verbose){cat(paste0('\nUpdating neighbours list with ', di, ' order neighbours...'))}
        if(di == 1){
          nns <- nns_raw
          di = di + 1
          next
        }
        nns <- lapply(nns, function(x){
          return(unique(unlist(nns_raw[x])))
        })
        di = di + 1
      }
      if(DEGREE != 1){
        nns <- lapply(1:length(nns), function(j){
          return(nns[[j]][(nns[[j]]!=j)])
        })
      }
      
      
      frequency_matrix <- do.call(rbind, lapply(1:length(nns), function(i){
        
        if(verbose){ 
          if(i %% 1000 == 0 | i == 1){ 
            cat(paste0(i, ' of ', length(nns), '...'))
          }
        }
        
        ## Removal of convex hull
        nnvector <- nns[[i]]
        if(i %in% convex_hull$i){
          nnvector <- nnvector[!(nnvector %in% convex_hull$i)]
        }
        
        ## Cutting by distance
        if(!is.null(delaunayTriangulationDistanceThreshold) & 
           !is.na(delaunayTriangulationDistanceThreshold) &
           !is.infinite(delaunayTriangulationDistanceThreshold)){
          ref_coord <- all_coords[i]
          coords <- all_coords[nnvector]
          dists <- Mod(coords - ref_coord)
          nnvector <- nnvector[dists < delaunayTriangulationDistanceThreshold]
        }
        
        values <- tabulate( CLB[nnvector] )
        if(length(values) < max(CLB)){
          values <- c(values, rep(0, max(CLB) - length(values)))
        }
        
        return(values)
      }))
      
      colnames(frequency_matrix) <- 
        paste0('DT', sprintf(paste0('%0', max(nchar(delaunayNNDegrees)), 'd'), DEGREE), '_', levels(factor(cell_label)))
      DTS[[j]] <- frequency_matrix 
    }
    DTS <- do.call(cbind, DTS)
    OUT[['DT']] <- DTS
  }
  ###
  
  
  ###
  if(!is.null(euclideanDistances)){
    
    euDists <- sort(unique(euclideanDistances))
    if(verbose){ cat('\n\nGetting nearest neighbours within specific euclidean distances...') }
    
    FM <- matrix(0, nrow = length(all_coords), ncol = max(CLB) * length(euDists))
    
    for(i in 1:length(all_coords)){
      if(verbose){
        if(i %% 1000 == 0 | i == 1){
          cat(paste0(i, ' of ', length(all_coords), '...'))
        }
      }
      dists <- Mod(all_coords - all_coords[i])
      reflab <- CLB[i] + c( (1:length(euDists)) - 1) * max(CLB)
      for(j in 1:length(reflab)){
        FM[,reflab[j]] <- FM[,reflab[j]] + (dists < euDists[j])
      }
    }
    
    colns <- sapply(1:length(euDists), function(j){
      paste0('EUD', 
             sprintf( paste0('%0', max(nchar(euclideanDistances)), 'd'), euDists[j]), 
             '_', levels(factor(cell_label)))
    })
    colnames(FM) <- as.vector(colns) 
    OUT[['FM']] <- FM
  }
  ###
  
  OUT <- do.call(cbind, OUT)
  return(OUT)
}