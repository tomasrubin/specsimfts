
##########################################################################################################
## integrate to get-h lag covariance operators
CKL_covlagh_operator <- function(harmonic_eigenvalues, harmonic_eigenfunctions, lag, n_grid, n_pc, n_grid_freq=1000){
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  
  integrand <- function(omega){ 
    ret <- array(0, dim=c(n_grid,n_grid))
    for (n in 1:n_pc){
      ret <- ret + Re( harmonic_eigenvalues(omega, n) * outer( harmonic_eigenfunctions(omega,n,grid), harmonic_eigenfunctions(omega,n,grid)) * exp(1i*omega*lag) )
    }
    return(ret)
  }
  
  # numeric integration
  int <- array(0, dim=c(n_grid,n_grid))
  for (k in 1:n_grid_freq){
    omega <- pi*k/n_grid_freq # 0..pi
    int <- int + integrand(omega) *pi /n_grid_freq
  }
  
  
  return( Re(int + t(int)) )
  
  
}