##########################################################################################################
## simulation in Cramer-Karhunen-Loeve representation

spec_density_simulate <- function(spec_density, t_max, n_grid, n_pc, seed_number = NULL, include_freq_zero = F){
  
  ## random seed if assigned
  if (!is.null(seed_number)){ set.seed(seed_number) }
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  
  # prepare the ts I'm generating
  t_half <- ceiling(t_max/2)
  t_max <- t_half*2
  ts <- c(1:t_half, t_max)
  
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
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # get Vs
  vs <- matrix(0, n_grid, t_max)
  #pb <- txtProgressBar(style = 3)
  for (ii in 1:(t_half+1)){
    # prepare eigenvalue and eigenfunction of Sigma
    omega <- (2*ii*pi)/t_max
    spec_density_svd <- svd(spec_density(omega, grid_matrix,t(grid_matrix)), nu=n_pc, nv=0)
    
    for (n in (1:n_pc)){
      
      # split by if we're at zero freq or not
      if (ii == t_half + 1){
        # this ii corresponds to t_max = zero frequency
        if (include_freq_zero){
          phi <- spec_density_svd$u[,n]
          vs[,t_max] <- vs[,t_max] + sqrt(spec_density_svd$d[n]) * phi * zs[n,ii]
        }
        # do nothing
      } else {
        # all the other ii, going from 1 to t_half
        
        # function
        phi <- spec_density_svd$u[,n]
        
        vs[,ii] <- vs[,ii] + sqrt(spec_density_svd$d[n]) * phi * zs[n,ii]
        
      }
    }
  }
  #close(pb)
  
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