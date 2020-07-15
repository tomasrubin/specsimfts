#' Simulate functional time series sample when the dynamics is defined as a filtered white noise through the white noise covariance operator and the filter is given by its frequency response function Theta.
#' 
#'
#' @title Simulate a FTS defined as a filter of a white noise
#' @param theta The frequency response operator Theta(omega) of the filter used for the definition of filtered white noise. A function of two variables, \code{omega} and \code{f}, where \code{f} is a vector (discretisation of a function) on which the operator Theta at frequency \code{omega} is applied onto. See the example bellow to inspect how you can define Theta(omega). The functions \code{\link{rank_one_tensor}} and \code{\link{kernel_operator}} might be useful for the definition of Theta. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi].
#' @param t_max Time horizon to be simulated. Must be an even number, otherwise it is increased by one.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1] where the FTS is to be simulated.
#' @param n_pc The number of eigenfunctions of sigma to be used for the simulation. Setting n_pc=n_grid is recommended as there is hardly any computational gain when n_pc is smaller. In case the sigma is defined as finite rank operator (through lists \code{sigma_eigenvalues} and \code{sigma_eigenfunctions}) and \code{n_pc} is higher than this finite rank, \code{n_pc} is automatically discresed to match this finite rank.
#' @param seed_number The random seed inicialization for the simulation. The value \code{NULL} means no inicialization
#' @param sigma The covariance operator of the white noise innovation. A function of two variables, \code{x} and \code{y}, returns the value of the covariance kernel evaluated at (\code{x},\code{y}).
#' @param sigma_eigenvalues Alternatively, you can define the white noise innovation covariance operator through its eigendecomposition, in which case supply the arguments \code{sigma_eigenvalues} and \code{sigma_eigenfunctions} but leave out \code{sigma} (or set it \code{NULL}). The eigendecomposition can be defined either as (i) a list of finite number of eigenvalues/eigenfunctions, or (ii) a function returning the eigenvalue/eigenfunction of given order. In the case (i), the \code{sigma_eigenvalues} parameter is simply a vector of eigenvalues. In the case (ii), the \code{sigma_eigenvalues} parameter is a function of the variable \code{n} which returns the \code{n}-th eigenvalue. See the example bellow for the implementation of the two cases.
#' @param sigma_eigenfunctions See \code{sigma_eigenvalues} before. In the case (i), the \code{sigma_eigenfunctions} parameter is a list of functions with the same length as the parameter \code{sigma_eigenvalues}. Each function is a function of the variable \code{x} that returns the value of the eigenfunction at point \code{x}. The order of the functions in the list determine the order of eigenfunctions. In the case (ii), the \code{sigma_eigenfunctions} parameter is a function of two variables, \code{n} and \code{x}, that returns the value of the \code{n}-th eigenfunction at the point \code{x}.
#' @param include_freq_zero If set \code{TRUE}, the zero frequency is included for simulation in the spectral domain. Set \code{FALSE} for processes with singularity at frequency zero, e.g. the long-range dependent FARFIMA(p,d,q) process.
#' @return functional time series sample, matrix of size (\code{n_grid},\code{t_max})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{filter_covlagh_operator}}, \code{\link{rank_one_tensor}}, \code{\link{kernel_operator}}
#' @examples
#' # define the white noise covariance operator (Brownian motion)
#' sigma <- function(x,y) { pmin(x,y) }
#'
#' # # Alternatively defined the sigma covariance through its eigendecomposition, and supply it to the function 'filter_simulate'
#' # sigma_eigenvalues <- function(n) { 1/((n-0.5)*pi)^2 }
#' # sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin((n-0.5)*pi*x) }
#' 
#' # # Alternatively, define the innovation covariance as low-rank (through a list)
#' # sigma_eigenvalues <- c(1, 0.6, 0.3, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)
#' # sigma_eigenfunctions <- list(
#' #   function(x){ sin(2*pi*x) },
#' #   function(x){ cos(2*pi*x) },
#' #   function(x){ sin(4*pi*x) },
#' #   function(x){ cos(4*pi*x) },
#' #   function(x){ sin(6*pi*x) },
#' #   function(x){ cos(6*pi*x) },
#' #   function(x){ sin(8*pi*x) },
#' #   function(x){ cos(8*pi*x) },
#' #   function(x){ sin(10*pi*x) },
#' #   function(x){ cos(10*pi*x) }
#' # )
#' 
#' # define filter
#' theta <- function(omega,f){
#' 2*f+
#'  1i*rev(f) +
#'  omega*cumsum(f)/length(f) +
#'  rank_one_tensor( function(x){sin(x)}, function(x){exp(x)}, f ) +
#'  kernel_operator( function(x,y){sin(omega+x+2*y)}, f )
#' }
#' 
#' # simulation setting
#' t_max <- 1000
#' n_grid <- 101
#' 
#' # simulate in the spectral domain
#' fts_x <- filter_simulate(theta, t_max, n_grid, sigma=sigma)
#' 
#' # # Alternatively simulate with the known eigendecomposition of sigma
#' # fts_x <- filter_simulate(theta, t_max, n_grid, sigma_eigenfunctions = sigma_eigenfunctions, sigma_eigenvalues=sigma_eigenvalues)
#' 
#' # plot the first curve
#' plot(fts_x[,1], type='l')
#' 
#' @export


filter_simulate <- function(theta, t_max, n_grid, n_pc=n_grid, seed_number = NULL, sigma=NULL, sigma_eigenvalues=NULL, sigma_eigenfunctions=NULL, include_zero_freq=F){
  
  ## random seed if assigned
  if (!is.null(seed_number)){ set.seed(seed_number) }
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # prepare the ts I'm generating
  t_half <- ceiling(t_max/2)
  t_max <- t_half*2
  ts <- c(1:t_half, t_max)
  
  # if provided analytic eigenfunctions or eigenvalues
  if (!is.null(sigma_eigenvalues) & !is.null(sigma_eigenfunctions)){
    
    # can be either function or list
    if (is.list(sigma_eigenvalues)){
      sigma_eigenfunctions_eval <- sigma_eigenfunctions
      for (n in 1:length(sigma_eigenfunctions)){
        sigma_eigenfunctions_eval[[n]] <- sigma_eigenfunctions(grid)
      }
      sigma_eig <- list(values=sigma_eigenvalues, vectors=sigma_eigenfunctions_eval)
    } else {
      # they're functions
      sigma_eigenvalues <- numeric(n_pc)
      sigma_eigenfunctions_eval <- matrix(NaN, nrow=n_grid, ncol=n_pc)
      for (n in 1:n_pc){
        sigma_eigenvalues[n] <- sigma_eigenvalues(n)
        sigma_eigenfunctions_eval[,n] <- sigma_eigenfunctions(n,grid)
      }
      sigma_eig <- list(values=sigma_eigenvalues, vectors=sigma_eigenfunctions_eval)
    }
    
  } else {
  # eigendecomposition of sigma, if not provided
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
    omega <- 0
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