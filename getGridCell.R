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
  
  df$Categorize <- paste0(df$gridGroup, "_", df$cell_type)
  # get a frequancy table of cell-types x grids of size g_n
  freq_table <- table(df$cell_type, df$gridGroup)
  # normalize so that each grid adds up to 1 
  normalized_data <- t(t(freq_table) / colSums(freq_table))
  # rename the grid
  colnames(normalized_data) <- paste0(seq(1, ncol(normalized_data)), "_grid")
  return(normalized_data)
}
