###################################################################################################
## simulation in a given filtration
filtration_simulate <- function(sigma, theta, t_max, n_grid, n_pc=n_grid, seed_number = NULL, sigma_eig = NaN, include_zero_freq=F){
  
  ## random seed if assigned
  if (!is.null(seed_number)){ set.seed(seed_number) }
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # prepare the ts I'm generating
  t_half <- ceiling(t_max/2)
  t_max <- t_half*2
  ts <- c(1:t_half, t_max)
  
  ## eigen decomposition of sigma, if not provided
  if (!is.list(sigma_eig)){ 
    sigma_grid <- sigma(grid_matrix, t(grid_matrix))
    sigma_svd <- svd( sigma_grid, nu=n_grid, nv = 0)
    sigma_eig <- list(values = sigma_svd$d, vectors = sigma_svd$u)
  }
  
  # adjust n_pc if higher then the rank of sigma
  if (n_pc > length(sigma_eig$values)){
    n_pc <- length(sigma_eig$values)
  }
  
  ## get variables Z
  
  # remember, t_half and t_max are real
  zs <- matrix(0, n_pc, t_half+1 )
  
  # generate zs, cycle through pc first. this is because of comparability across different n_pc with the same seed
  for (pc in 1:n_pc){
    # complex zs
    zs[pc,1:(t_half-1)] <- rnorm(t_half-1) + 1i * rnorm(t_half-1)
    # real zs
    zs[pc,t_half] <- 2*rnorm(1)
    zs[pc,t_half+1] <- 2*rnorm(1)
  }
  
  
  ## get variables V
  
  # get Vs
  vs <- matrix(0, n_grid, t_max)
  #pb <- txtProgressBar(style = 3)
  for (ii in 1:(t_half+1)){
    omega <- (2*ii*pi)/t_max
    
    # function to be plugged into Theta
    f <- numeric(n_grid)
    for (n in 1:n_pc){
      f <- f + sigma_eig$vectors[,n] * (sqrt(sigma_eig$values[n] / (2*pi)) * zs[n,ii])  # the constant (2*pi) here is due to the fact that the white noise has the spectral density (1/2*pi)*Sigma
    }
    # f <- sigma_eig$vectors %*% (sqrt(sigma_eig$values) * zs[,ii])
    
    # evaluate theta
    vs[,ii] <- theta(omega, f)
    
  }
  # zero freq
  if (include_zero_freq){
    f <- numeric(n_grid)
    for (n in 1:n_pc){
      f <- f + sigma_eig$vectors[,n] * (sqrt(sigma_eig$values[n]/ (2*pi)) * zs[n,t_half+1]) 
    }
    vs[,t_max] <-  theta(omega, f)
  }
  
  # mirror Vs t_half+1,...,t_max-1
  for (ii in 1:(t_half-1)){
    vs[,t_max-ii] <- Conj(vs[,ii])
  }
  
  ## iFFT to temporal domain
  
  # rearange the array so zero freq is at the beginning - that's what iFFT needs
  vs_for_ifft <- matrix(0, n_grid, t_max)
  vs_for_ifft[,2:t_max] <- vs[,1:(t_max-1)]
  vs_for_ifft[,1] <- vs[,t_max] # zero freq
  
  return( sqrt(pi/t_max) * Re(t(mvfft( t(vs_for_ifft), inverse = T ))) )
  
}