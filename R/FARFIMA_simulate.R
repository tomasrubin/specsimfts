###################################################################################################
## simulation in a given filtration
FARFIMA_simulate <- function(FARFIMA_pars, t_max, n_grid, seed_number=NULL, hybrid_ar=T, burnin=100){
  
  ## random seed if assigned
  if (!is.null(seed_number)){ set.seed(seed_number) }
  
  # evaluate operators on grid
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # check existence of sigma
  if(!is.null(FARFIMA_pars[["sigma"]])){
    # sigma is defined
    if (is.character(FARFIMA_pars$sigma)){
      if (FARFIMA_pars$sigma == "bm"){ # if sigma is given by string "bm", it means it's the brownian motion
        go_check_eigenvalues <- F
        FARFIMA_pars$sigma <- function(x,y){pmin(x,y)}
        sigma_eig <- BM_eig(n_grid)
      } else {
        # if it's string but not "bm" => error
        stop("Error: FARFIMA_pars$sigma is string but not 'bm' which is the only acceptable string value") 
      }
    } else {
      # if sigma is defined but not "bm", go first check the definition of eigenvalues, if they're defined, use them instead of numerical SVD
      go_check_eigenvalues <- T
    }
  } else {
    # sigma is NOT defined
    go_check_eigenvalues <- T
  }
  
  if (go_check_eigenvalues){
    # if assigned eigen functions of sigma
    if (!is.null(FARFIMA_pars[["sigma_eigenfunctions"]])){
      # evaluate eigenfunctions
      sigma_eigenfunctions_eval = matrix(0, ncol=length(FARFIMA_pars$sigma_eigenfunctions), nrow=n_grid)
      for (ii in 1:length(FARFIMA_pars$sigma_eigenfunctions)){
        sigma_eigenfunctions_eval[,ii] = FARFIMA_pars$sigma_eigenfunctions[[ii]](grid)
      }
      
      # save into list
      sigma_eig <- list(
        values = FARFIMA_pars$sigma_eigenvalues,
        vectors = sigma_eigenfunctions_eval
      )
      
      # if not specified, create the function for sigma
      if(!is.null(FARFIMA_pars[["sigma"]])){
        FARFIMA_pars$sigma <- kernel_from_eig( FARFIMA_pars$sigma_eigenvalues, FARFIMA_pars$sigma_eigenfunctions )
      }
    } else {
      # numerically evaluate the eigendecomposition of sigma
      sigma_svd <- svd( FARFIMA_pars$sigma( grid_matrix,t(grid_matrix) ), nu = n_grid, nv = 0 )
      
      # save into list
      sigma_eig <- list(
        values = sigma_svd$d,
        vectors = sigma_svd$u
      )
    }
  }
  
  n_pc <- length(sigma_eig$values)
  
  # save model order
  ar_order <- length(FARFIMA_pars$operators_ar)
  ma_order <- length(FARFIMA_pars$operators_ma)
  
  
  ### evaluate the MA and AR operators
  operators_ar_eval <- FARFIMA_pars$operators_ar
  operators_ma_eval <- FARFIMA_pars$operators_ma
  
  if (ar_order>0){
    for (j in 1:ar_order){
      operators_ar_eval[[j]] <- FARFIMA_pars$operators_ar[[j]](grid_matrix,t(grid_matrix))
    }
  }
  if (ma_order>0){
    for (j in 1:ma_order){
      operators_ma_eval[[j]] <- FARFIMA_pars$operators_ma[[j]](grid_matrix,t(grid_matrix))
    }
  }
  
  # check if I'm doing burnin and then increase t_max by the burnin
  if ((ar_order>0) & (hybrid_ar)){
    t_max <- t_max + burnin
  } else {
    hybrid_ar <- FALSE
  }
  
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
  
  # get Vs
  vs <- matrix(0, n_grid, t_max)
  
  for (ii in 1:(t_half+1)){
    omega <- (2*ii*pi)/t_max
    
    # function to be plugged into Theta
    f <- numeric(n_grid)
    for (n in 1:n_pc){
      f <- f + sigma_eig$vectors[,n] * (sqrt(sigma_eig$values[n] / (2*pi)) * zs[n,ii])  # the constant (2*pi) here is due to the fact that the white noise has the spectral density (1/2*pi)*Sigma
    }
    # f <- sigma_eig$vectors %*% (sqrt(sigma_eig$values) * zs[,ii])
    
    # evaluate theta
    # apply MA part (if exsits)
    if (ma_order>0){
      ma_part <- diag(n_grid)
      for (j in 1:ma_order){
        ma_part <- ma_part + operators_ma_eval[[j]] * exp(-1i*omega*j) / n_grid
      }
      f <- ma_part %*% f
    }
    
    # apply AR part (if exists)
    if ((ar_order>0) & (!hybrid_ar)){ # only if I'm NOT doing the AR in the time domain
      ar_part <- diag(n_grid)
      for (j in 1:ar_order){
        ar_part <- ar_part + operators_ar_eval[[j]] * exp(-1i*omega*j) / n_grid
      }
      f <- solve(ar_part, f)
    }
    
    # apply fractional integration
    vs[,ii] <- f * ( 2 * sin(omega/2) )^(-FARFIMA_pars$fractional_d)
    
    
  }
  # zero freq
  if (FARFIMA_pars$fractional_d <= 0){
    f <- numeric(n_grid)
    for (n in 1:n_pc){
      f <- f + sigma_eig$vectors[,n] * (sqrt(sigma_eig$values[n]/ (2*pi)) * zs[n,t_half+1]) 
    }
    
    # evaluate theta
    # apply MA part (if exsits)
    if (ma_order>0){
      ma_part <- diag(n_grid)
      for (j in 1:ma_order){
        ma_part <- ma_part + operators_ma_eval[[j]] * exp(-1i*omega*j) / n_grid
      }
      f <- ma_part %*% f
    }
    
    # apply AR part (if exists)
    if ((ar_order>0) & (!hybrid_ar)){ # only if I'm NOT doing the AR in the time domain
      ar_part <- diag(n_grid)
      for (j in 1:ar_order){
        ar_part <- ar_part + operators_ar_eval[[j]] * exp(-1i*omega*j) / n_grid
      }
      f <- solve(ar_part, f)
    }
    
    # apply fractional integration
    vs[,t_max] <-  f * ( 2 * sin(omega/2) )^(-FARFIMA_pars$fractional_d)
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
  
  fts_x <- sqrt(pi/t_max) * Re(t(mvfft( t(vs_for_ifft), inverse = T )))
  
  if (hybrid_ar){
    # apply the AR recursion in the time domain
    fts_x <- apply_AR_part(fts_x, FARFIMA_pars$operators_ar)
    fts_x <- fts_x[,(burnin+1):ncol(fts_x)] 
    return( fts_x )
  } else {
    # just simply return fts_x
    return( fts_x )
  }
  
}