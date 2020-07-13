#' Simulate the FARFIMA(p,d,q) process, the functional autoregressive fractionally integrated moving average process.
#' 
#'
#' @title Simulate the FARFIMA(p,d,q) process in the spectral domain
#' @param FARFIMA_pars The list of the parameters for the FARFIMA(p,d,q) process. Must contain fields: (i) \code{fractional_d}, a real number in the open interval (-0.5,0.5) controling the fractional integration degree. \code{fractional_d} being positive corresponds to long-rande dependence behaviour. (ii) \code{operators_ar}, the list of length 'p' the order of the autoregressive part. The autoregressive operators are considered to be integral operators defined through their kernels which are saved as the elements of the list \code{operators_ar} as functions of two variables, \code{x} and \code{y}, returning the value of the kernel at point (\code{x},\code{y}). In case of degenerate autoregressive part define \code{operators_ar} as an empty list. (iii) \code{operators_ma}, the list of length 'q', the order of the moving average part. Just like \code{operators_ar} its a liks of functions - the kernels of the moving average operators. (iv) The covariance opperator of the stochastic innovation process can be defined either through (iv-a) its kernel,  (iv-b) finite rank eigendecomposition, (iv-c) infinite rank decomposition. In the case (iv-a), define \code{sigma} as a function of two variables \code{x} and \code{y}, returning the value of the covariance kernel at point (\code{x},\code{y}). In the case (iv-b), define the elements \code{sigma_eigenvalues} as a vector of finitely many eigenvalues and \code{sigma_eigenfunctions} as a list of the same length as \code{sigma_eigenvalues} with each element being a function of variable \code{x} returning the value of that eigenfunction at point \code{x}. In the case (iv-c), define the elements \code{sigma_eigenvalues} as a function of the variable \code{n} returning the \code{n}-th eigenvalue and the element \code{sigma_eigenfunctions} as a function of two variables, \code{n} and \code{x}, returning the value of the \code{n}-th eigenfunctions at point \code{x}. See the example bellow for some examples on how to set up \code{FARFIMA_pars}.
#' @param t_max Time horizon to be simulated. Must be an even number, otherwise it is increased by one.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1] where the FTS is to be simulated.
#' @param seed_number The random seed inicialization for the simulation. The value \code{NULL} means no inicialization
#' @param hybrid_ar If set \code{TRUE} (default), the method first simulates the corresponding FARFIMA(0,d,q) process and then applies the autoregressive filter in the temporal domain. Such runtime is much faster than the fully spectral method (\code{hybrid_ar = FALSE}) in which case the method needs to solve a system of linear equation at each frequency.
#' @param burnin If hybrid_ar=TRUE set how long is the burn-in period for the autoregressive part (100 by default) in the temporal domain.
#' @return functional time series sample, matrix of size (\code{n_grid},\code{t_max})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{FARFIMA_covlagh_operator}} \code{\link{FARFIMA_test_stationarity}}
#' @examples
#' # (i) fractional integration
#' fractional_d <- 0 # in the open interval (-0.5, 0.5), positive number means long-range dependence
#' 
#' # (ii) autoregressive operators
#' operators_ar <- list(
#' function(x,y){ 0.3*sin(x-y) },
#' function(x,y){ 0.3*cos(x-y) },
#' function(x,y){ 0.3*sin(2*x) },
#' function(x,y){ 0.3*cos(y) }
#' )
#' # operators_ar <- list() # use empty list for degenerate AR part
#' 
#' # (iii) moving average kernels
#' # you can put here arbitrary long list of operators
#' operators_ma <- list(
#' function(x,y){ x+y },
#' function(x,y){ x },
#' function(x,y){ y }
#' )
#' # operators_ma <- list() # use empty list for degenerate MA part
#' 
#' # (iv-b) covariance of the inovation defined through eigenvalues and eigenfunctions
#' # you can put here arbitrary long lists but their lenghts should match
#' sigma_eigenvalues <- c(1, 0.6, 0.3, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)
#' sigma_eigenfunctions <- list(
#' function(x){ sin(2*pi*x) },
#' function(x){ cos(2*pi*x) },
#' function(x){ sin(4*pi*x) },
#' function(x){ cos(4*pi*x) },
#' function(x){ sin(6*pi*x) },
#' function(x){ cos(6*pi*x) },
#' function(x){ sin(8*pi*x) },
#' function(x){ cos(8*pi*x) },
#' function(x){ sin(10*pi*x) },
#' function(x){ cos(10*pi*x) }
#' )
#' 
#' # # (iv-c) innovation covariance operator (Brownian motion)
#' # sigma_eigenvalues <- function(n) { 1/((n-0.5)*pi)^2 }
#' # sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin((n-0.5)*pi*x) }
#' 
#' # put the parameters into one list
#' FARFIMA_pars <- list(fractional_d=fractional_d, operators_ar=operators_ar, operators_ma=operators_ma, sigma_eigenvalues=sigma_eigenvalues,sigma_eigenfunctions=sigma_eigenfunctions)
#' 
#' # # (iv-a) Alternatively, define the kernel of the white noise innovation.
#' # sigma <- function(x,y) { pmin(x,y) } # Brownian motion
#' # FARFIMA_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma, sigma=sigma)
#' 
#' 
#' # simulate trajectory
#' if (FARFIMA_test_stationarity(FARFIMA_pars)){
#'  # fully spectral approach. if AR part is non-degenerate, the simulation involves solving a system of liear equations at each frequency
#'  fts_x <- FARFIMA_simulate(FARFIMA_pars, t_max, n_grid, hybrid_ar = F)
#'    
#'  # # hybrid simulation method
#'  # fts_x <- FARFIMA_simulate(FARFIMA_pars, t_max, n_grid, hybrid_ar = T)
#'    
#'  # display the first curve
#'  plot(fts_x[,1], type='l')
#' }
#' 
#' @export

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
      
      if (is.list(FARFIMA_pars[["sigma_eigenfunctions"]])){
        # low-rank specification
        
        # evaluate eigenfunctions
        sigma_eigenfunctions_eval <- matrix(0, ncol=length(FARFIMA_pars$sigma_eigenfunctions), nrow=n_grid)
        for (ii in 1:length(FARFIMA_pars$sigma_eigenfunctions)){
          sigma_eigenfunctions_eval[,ii] <- FARFIMA_pars$sigma_eigenfunctions[[ii]](grid)
        }
        
        # save eigenvalues
        sigma_eigenvalues_eval <- FARFIMA_pars$sigma_eigenvalues
      } else {
        # functions for eigenvalues, eigenfunctions
        
        # evaluate eigenfunctions
        n_pc <- n_grid
        sigma_eigenfunctions_eval <- matrix(0, ncol=n_pc, nrow=n_grid)
        sigma_eigenvalues_eval <- numeric(n_pc)
        for (ii in 1:n_pc){
          sigma_eigenfunctions_eval[,ii] <- FARFIMA_pars$sigma_eigenfunctions(ii,grid)
          sigma_eigenvalues_eval[ii] <- FARFIMA_pars$sigma_eigenvalues(ii)
        }
      }
      
      # save into list
      sigma_eig <- list(
        values = sigma_eigenvalues_eval,
        vectors = sigma_eigenfunctions_eval
      )
      
      
      # if not specified, create the function for sigma
      if(is.null(FARFIMA_pars[["sigma"]])){
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
        ar_part <- ar_part - operators_ar_eval[[j]] * exp(-1i*omega*j) / n_grid
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
        ar_part <- ar_part - operators_ar_eval[[j]] * exp(-1i*omega*j) / n_grid
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