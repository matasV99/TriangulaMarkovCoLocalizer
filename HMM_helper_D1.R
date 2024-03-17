# helper functions for HMM testing 
easypackages::libraries("tidyverse", "data.table", "future", "fasano.franceschini.test", "scales", "mHMMbayes")
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
  idx <- x + (x) * (y-min(y))
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
      reflab <- c(CLB[i], CLB[i] + max(CLB))
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


getGridCell <- function(
    ## INPUT
  x, ## X coordinates
  y, ## Y coordinates
  cell_label, ## Cell type labels
  g_n = NULL # length of grid edge
){
  
  if(length(x) != length(y) | length(x) != length(cell_label) ){
    cat('\n\nERROR: X, Y, AND CELL_LABELS REQUIRE THE SAME NUMBER OF ENTRIES!')
    
  }
  
  if(is.null(g_n)){
    cat("\n\nERROR: Specify distance length of grid edge!")
  }
  
  df <- data.frame(x = x, y = y)
  df$cell_type <- cell_label
  
  df$gridGroup <- interaction(cut(df$x, breaks=seq(0, max(df$x) + g_n, by=g_n)),
                              cut(df$y, breaks=seq(0, max(df$y) + g_n, by=g_n)), sep=";")
  return(df)
}


getCelltypeProbs <- function(df, celltype_vec = celltype_vec) {

  dt_cols <- names(df)[grepl("DT", names(df))]
  df_mean <- df %>%
    group_by(cell_type) %>%
    summarise_at(dt_cols, mean) %>%
  mutate(across(where(is.numeric), ~./rowSums(across(where(is.numeric))))) %>%
  mutate(across(everything(), ~replace(., is.nan(.), 0)))

  # pad by cell_types with no representation in grid 
  missing_celltype <- celltype_vec[!(celltype_vec %in% df_mean$cell_type)]

  df_export <- df_mean %>% bind_rows(
    tibble(cell_type = missing_celltype)) %>%
    mutate(across(everything(), ~replace(., is.na(.), 0)))

  return(df_export)
}

calc_ks_dist4 <- function(grid_list) {
  num_grids <- length(grid_list)
  
  dist_mat <- matrix(NA, nrow = num_grids, ncol = num_grids)
  
  # Create a list of all combinations of i and j
  ij_combinations <- combn(1:num_grids, 2, simplify = FALSE)
  
  # Define a function to calculate the distance for a given pair of i and j
  calc_distance <- function(ij) {
    i <- ij[1]
    j <- ij[2]
    
    output <- fasano.franceschini.test(as.matrix(grid_list[[i]][-c(1:5)]), 
                                       as.matrix(grid_list[[j]][-c(1:5)]),
                                       nPermute=0 ,seed=100, threads=35, verbose=TRUE)
    dist <- unname(output[[1]])
    
    list(i = i, j = j, dist = dist)
  }
  
  # Use mclapply to apply the function to each pair of i and j in parallel
  results <- mclapply(ij_combinations, calc_distance, mc.cores = detectCores())
  
  # Fill in the distance matrix with the results
  for(result in results) {
    dist_mat[result$i, result$j] <- result$dist
    dist_mat[result$j, result$i] <- result$dist
  }
  
  return(dist_mat)
}

calc_px_dist <- function(grid_list) {
  num_grids <- length(grid_list)
  
  # Initialize a matrix to store the centroids
  centroid_mat <- matrix(NA, nrow = num_grids, ncol = 2)
  
  # Calculate the mean x and y coordinates for each grid
  for(i in 1:num_grids) {
    centroid_mat[i, 1] <- mean(grid_list[[i]]$x)
    centroid_mat[i, 2] <- mean(grid_list[[i]]$y)
  }
  
  # Calculate the Euclidean distance between the centroids
  dist_mat <- as.matrix(dist(centroid_mat, method = 'euclidean'))
  
  return(dist_mat)
}
