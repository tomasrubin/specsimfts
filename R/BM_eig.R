###################################################################################################
## analytic SVD decomposition od the Brownian motion kernel
BM_eig <- function(n_grid){
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  
  sigma_eig <- list()
  sigma_eig$vectors <- matrix(0,nrow=n_grid,ncol=n_grid)
  sigma_eig$values <- numeric(n_grid)
  for (n in 1:n_grid){
    sigma_eig$vectors[,n] <- sqrt(2) * sin( (n-0.5)*pi*grid )
    sigma_eig$values[n]  <- 1/(((n-0.5)*pi)^2) 
  }
  
  return(sigma_eig)
}
