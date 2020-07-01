
##########################################################################################################
## integrate to get-h lag covariance operators
spec_density_covlagh_operator <- function(spec_density, lag, n_grid, n_grid_freq=2000){
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  int <- array(0, dim=c(n_grid,n_grid))
  for (k in 1:n_grid_freq){
    omega <- pi*k/n_grid_freq # 0..pi
    int <- int + spec_density(omega,grid_matrix,t(grid_matrix))* exp(1i*omega*lag) *pi /n_grid_freq 
  }
  
  return( Re(int + t(int)) )
  
}