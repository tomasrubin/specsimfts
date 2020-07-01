
##########################################################################################################
## integrate to get-h lag covariance operators
filtration_covlagh_operator <- function(sigma, theta, lag, n_grid, n_grid_freq=2000){
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  int <- array(0, dim=c(n_grid,n_grid))
  for (k in 1:n_grid_freq){
    omega <- pi*k/n_grid_freq # 0..pi
    
    # construct theta discretization as a kernel
    theta_eval <- matrix(0,nrow=n_grid,ncol=n_grid)
    for (ii in 1:n_grid){
      theta_eval[,ii] <- theta(omega, replace(numeric(n_grid), ii, 1))
    }
    
    spec_density <- 1/(2*pi) * theta_eval %*% sigma(grid_matrix,t(grid_matrix)) %*% Conj(t(theta_eval))
    int <- int + spec_density * exp(1i*omega*lag) *pi /n_grid_freq 
  }
  
  return( Re(int + t(int)) )
  
}