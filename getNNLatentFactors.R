getNNLatentFactors <- function(
    
  ## INPUT
  NNMatrix,
  method = 'nmf',
  
  ## Global params
  k = 5,
  seed = 12345,
  verbose = TRUE,
  
  ## RcppML nmf() params
  tol = 1e-04,
  maxit = 200,
  L1 = c(0, 0),
  mask_zeros = FALSE,
  diag = TRUE,
  nonneg = TRUE,
  
  ## CorEx params
  eps=1e-5, 
  count='binarize',
  tree=TRUE
  
){
  
  ###
  seed = as.integer(seed)
  k = as.integer(k)
  maxit = as.integer(maxit)
  ###
  
  ###
  supported_methods <- c('nmf', 'corex', 'hmm')
  method = tolower(method)
  if(!(method %in% supported_methods)){
    cat("\n\nERROR: METHOD NOT SUPPORTED! ONLY AVAILABLE METHODS:")
    for(i in 1:length(supported_methods)){cat(paste0('\n', supported_methods[i]))}
    break
  }
  ###
  
  ###
  OUT <- list(
    'method' = method
  )
  ###
  
  ###
  if(method=='nmf'){
    if(verbose){ cat('\n\nRunning NMF approach...') }
    require(RcppML)
    set.seed(seed)
    nmfd <- RcppML::nmf(
      NNMatrix, 
      k=k,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      L1 = L1,
      seed = seed,
      mask_zeros = mask_zeros,
      diag = diag,
      nonneg = nonneg)
    cat('Getting cell-wise loadings...')
    w <- nmfd$w
    cat('Getting feature-wise loadings...')
    h <- nmfd$h
    colnames(w) <- paste0('NMF', sprintf(paste0('%0', nchar(k), 'd'), 1:k))
    colnames(h) <- colnames(nndf)
    OUT[['CELL_LOADINGS']] <- w
    OUT[['FEATURE_LOADINGS']] <- h
  }
  ###
  
  ###
  if(method=='corex'){
    require(reticulate)
    if(!reticulate::py_available()){
      cat('\n\nERROR: NO PYTHON ENVIRONMENT HAS BEEN LOADED! TRY:\n')
      cat( '\nlibrary(reticulate)' )
      cat( '\nuse_python( {PYTHON LOCATION} )' )
      cat( '\nuse_condaenv( {CONDA ENV NAME}, required = TRUE )')
      break
    }
    
    if(!py_module_available('corextopic.corextopic')){
      cat('\n\nERROR: IT SEEMS YOUR PYTHON ENVIRONMENT DOES NOT HAVE THE COREX MODULE! TRY:\n')
      cat('\nconda activate {CONDA ENV NAME}')
      cat('\npip install corextopic')
      break
    }
    
    if(verbose){ cat('\n\nRunning CorEx model...')}
    corextopic <- import('corextopic.corextopic')
    norm <- NNMatrix>0
    # cat('Building model..')
    model <- corextopic$Corex(
      n_hidden = k, 
      seed = seed,
      eps = eps, 
      count = count,
      tree= tree
    )

    # cat('Fitting model...')
    model$fit(  
      norm, words = colnames(NNMatrix), docs = 1:nrow(norm)
    )
    cat('Getting cell-wise loadings...')
    w <- list()
    for(i in 1L:k){
      # cat(paste0(i, ' of ', k, '...'))
      vals <- unlist(model$get_top_docs(n_docs = nrow(norm), topic = i-1L))
      idx <- as.integer(vals[c(T, F)])
      logprob <- as.numeric(vals[c(F, T)])
      logprob <- logprob[order(idx)]
      w[[i]] <- logprob
      rm(idx, vals, logprob)
    }
    w <- do.call(cbind, w)
    colnames(w) <- paste0('CX', sprintf(paste0('%0', nchar(k), 'd'), 1:k))
    
    cat('Getting feature-wise loadings...')
    lfs <- model$get_topics(n_words=ncol(NNMatrix))
    facs <- sapply(1:length(lfs), function(x){
      if(length(lfs[[x]]) == 0){
        return(NA)
      }
      val <- unlist(lfs[[x]])
      if(val[3]!=1){
        nval <- (rev(val)[c(F, F, T)])
      }else{
        nval <- (val[c(T, F, F)])
      }
      return(sort(nval))
    })
    h <- matrix(0, nrow = k, ncol = ncol(NNMatrix))
    colnames(h) <- colnames(NNMatrix)
    for(i in 1:k){
      h[i,facs[[i]]] <- 1
    }
    OUT[['CELL_LOADINGS']] <- w
    OUT[['FEATURE_LOADINGS']] <- h
  }
  ###
  
  ###
  if( method == 'hmm'){
    cat('ERROR: METHOD CURRENTLY UNSUPPORTED!')
    break
  }
  ###
  
  return(OUT)
}
