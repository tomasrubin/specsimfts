############################################################################################################################################
# test stationarity of FARFIMA
FARFIMA_test_stationarity <- function(FARFIMA_pars, n_grid){
  
  # evaluate operators on grid
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # save model order
  ar_order <- length(FARFIMA_pars$operators_ar)
  
  
  if (ar_order == 0){# if no AR part, always stationary
    return(TRUE)
  } else {
    
    # express AR operators
    operators_ar_eval <- FARFIMA_pars$operators_ar
    
    for (j in 1:ar_order){
      operators_ar_eval[[j]] <- FARFIMA_pars$operators_ar[[j]](grid_matrix,t(grid_matrix)) / n_grid
    }
    
    # compose the AR operator in the state-space
    composed_a <- matrix(0, nrow=ar_order*n_grid, ncol=ar_order*n_grid)
    for (ii in 1:ar_order){
      composed_a[ 1:n_grid , ((ii-1)*n_grid+1):(ii*n_grid)] <- operators_ar_eval[[ii]]
    }
    if (ar_order>1){
      for (ii in 1:(ar_order-1)){
        composed_a[ (ii*n_grid+1):((ii+1)*n_grid), ((ii-1)*n_grid+1):(ii*n_grid) ] <- diag(n_grid)
      }
    }
    
    m <- composed_a %^% 1000 # perform matrix power
    
    if (is.na(max(m))){
      return(F)
    } else {
      if (is.finite(m[1,1])){
        return( norm( m, type="2" ) < 1 )
      } else {
        return(F)
      }
    }
  }
}