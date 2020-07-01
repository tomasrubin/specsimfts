###################################################################################################
## simulation in a given filtration
apply_AR_part <- function(fts_x, operators_ar){
  
  # first "ar_order" elements of the ouput FTS are zeros
  
  # grid for evaluation
  n_grid <- nrow(fts_x)
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # shortcut for AR oder
  ar_order <- length(operators_ar)
  
  # do the AR recursion only if nondegenerate AR part
  if (ar_order > 0){
    
    t_max <- ncol(fts_x)
    fts_x_new <- matrix(0, nrow=nrow(fts_x), ncol=t_max)
    
    # express AR operators
    ar_evaluation <- operators_ar
    for (j in 1:ar_order){
      ar_evaluation[[j]] <- operators_ar[[j]](grid_matrix,t(grid_matrix)) 
    }
    
    # AR cycle
    for (t in (1+ar_order):t_max){
      fts_x_new[,t] <- fts_x[,t] # noise contribution
      for (j in 1:ar_order){
        fts_x_new[,t] <- fts_x_new[,t] - ar_evaluation[[j]] %*% fts_x_new[,t-j] / n_grid
      }
    }
    
    return(fts_x_new)
    
  } else {
    return(fts_x)
  }
  
}