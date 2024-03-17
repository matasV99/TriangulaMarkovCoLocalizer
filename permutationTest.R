permutationTest <- function(
    
  ## INPUT
  NNM,
  cell_labels,
  
  alternative = 'greater',
  nPerms = 1E3,
  seed = 12345,
  verbose = TRUE
  
){
  
  alternative = tolower(alternative)
  supported_alts <- c('greater', 'lower', 'two-tailed')
  if(!(alternative %in% supported_alts)){
    cat('\n\nERROR: UNSUPPORTED ALTERNATIVE! ONLY ACCEPTING:')
    for(i in 1:length(supported_alts)){ cat(paste0('\n', alternative, '...'))}
    break
  }
  
  if(is.null(dim(NNM))){
    NNMatrixVector <- NNM
    obs <- by(NNMatrixVector, cell_labels, sum)
    obs <- setNames(as.numeric(obs), names(obs))
    
    set.seed(seed)
    nulldist <- list()
    for(px in 1:nPerms){
      if(verbose){
        if(px==1 | px %% 100 == 0){
          cat(paste0(px, ' of ', nPerms, '...'))
        }
      }
      null_labs <- sample(cell_labels, size = length(cell_labels), replace = F)
      expx <- by(NNMatrixVector, null_labs, sum)
      expx <- setNames(as.numeric(expx), names(expx))
      nulldist[[px]] <- expx
    }
    nulldist <- do.call(cbind, nulldist)
    mean_exp <- apply(nulldist, 1, mean)
    sd_exp <- apply(nulldist, 1, sd)
    
    pvals <- sapply(1:nrow(nulldist), function(j){
      pval = 1 - sum( obs[j] < nulldist[j,]) / nPerms
      if(alternative == 'greater'){
        pval = 1 - pval
      }
      if(alternative == 'two-tailed'){
        pval = pmin(1 - pval, pval)
        pval = pval * 2
      }
      return(pval)
    })
    
    OUT <-  data.frame(
      'OBS' = obs,
      'OBS PERCENT' = round(100 * obs / sum(obs), digits=2), 
      'MEAN EXP' = mean_exp,
      'SD EXP' = sd_exp,
      'Z' = (obs - mean_exp)/sd_exp,
      'PVAL' = pvals,
      'SIG' = c('***', '**', '*', '.', '')[as.integer(cut(pvals, c(1, 0.1, 0.05, 0.01, 0.001, -1)))]
    )
    return(OUT)
  }
  
  start_time <- Sys.time()
  l_o <- list()
  for (c in unique(cell_labels)){
    NNMatrixVector <- NNM[,c]
    obs <- by(NNMatrixVector, cell_labels, sum)
    obs <- setNames(as.numeric(obs), names(obs))
    
    set.seed(seed)
    nulldist <- list()
    for(px in 1:nPerms){
      if(verbose){
        if(px==1 | px %% 100 == 0){
          cat(paste0(px, ' of ', nPerms, '...'))
        }
      }
      null_labs <- sample(cell_labels, size = length(cell_labels), replace = F)
      expx <- by(NNMatrixVector, null_labs, sum)
      expx <- setNames(as.numeric(expx), names(expx))
      nulldist[[px]] <- expx
    }
    nulldist <- do.call(cbind, nulldist)
    mean_exp <- apply(nulldist, 1, mean)
    sd_exp <- apply(nulldist, 1, sd)
    
    pvals <- sapply(1:nrow(nulldist), function(j){
      pval = 1 - sum( obs[j] < nulldist[j,]) / nPerms
      if(alternative == 'greater'){
        pval = 1 - pval
      }
      if(alternative == 'two-tailed'){
        pval = pmin(1 - pval, pval)
        pval = pval * 2
      }
      return(pval)
    })
    
    OUT <-  data.frame(
      'OBS' = obs,
      'OBS PERCENT' = round(100 * obs / sum(obs), digits=2), 
      'MEAN EXP' = mean_exp,
      'SD EXP' = sd_exp,
      'Z' = (obs - mean_exp)/sd_exp,
      'PVAL' = pvals,
      'SIG' = c('***', '**', '*', '.', '')[as.integer(cut(pvals, c(1, 0.1, 0.05, 0.01, 0.001, -1)))]
    )
    l_o[[c]] <- OUT
  }
  combined_df <- do.call("rbind", l_o)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Total time elapsed:", time_elapsed))
  return(combined_df)
}
