#' Numerically calculate the lag-h covariance operators for functional time series dynamics defined as a filter of a white noise. The calculation is done by numerically integrating the inverse formula, i.e. the spectral density multiplied by exp(-1i*lag*omega)
#'
#' @title Calculate the lag-h autocovariance operator for FTS defined as a filter of a white noise
#' @param theta The frequency response operator Theta(omega) of the filter used for the definition of filtered white noise. A function of two variables, \code{omega} and \code{f}, where \code{f} is a vector (discretisation of a function) on which the operator Theta at frequency \code{omega} is applied onto. See the example bellow to inspect how you can define Theta(omega). The functions \code{\link{rank_one_tensor}} and \code{\link{kernel_operator}} might be useful for the definition of Theta. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi].
#' @param lag The lag of the autocovariance to evaluate.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1]^2 for the operator kernel to evaluate.
#' @param sigma The covariance operator of the white noise innovation. A function of two variables, \code{x} and \code{y}, returns the value of the covariance kernel evaluated at (\code{x},\code{y}).
#' @param sigma_eigenvalues Alternatively, you can define the white noise innovation covariance operator through its eigendecomposition, in which case supply the arguments \code{sigma_eigenvalues} and \code{sigma_eigenfunctions} but leave out \code{sigma} (or set it \code{NULL}). The eigendecomposition can be defined either as (i) a list of finite number of eigenvalues/eigenfunctions, or (ii) a function returning the eigenvalue/eigenfunction of given order. In the case (i), the \code{sigma_eigenvalues} parameter is simply a vector of eigenvalues. In the case (ii), the \code{sigma_eigenvalues} parameter is a function of the variable \code{n} which returns the \code{n}-th eigenvalue. See the example bellow for the implementation of the two cases.
#' @param sigma_eigenfunctions See \code{sigma_eigenvalues} before. In the case (i), the \code{sigma_eigenfunctions} parameter is a list of functions with the same length as the parameter \code{sigma_eigenvalues}. Each function is a function of the variable \code{x} that returns the value of the eigenfunction at point \code{x}. The order of the functions in the list determine the order of eigenfunctions. In the case (ii), the \code{sigma_eigenfunctions} parameter is a function of two variables, \code{n} and \code{x}, that returns the value of the \code{n}-th eigenfunction at the point \code{x}.
#' @param n_grid_freq The grid points for the spectral density to evaluate at. Partition of [0,pi].
#' @return lag-h autocovariance operator, matrix of size (\code{n_grid},\code{n_grid})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{filter_simulate}}
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
#' # evaluation setting
#' lag <- 1 # change here to evaluate different lag-h autocovariance operator. put "lag <- 0" for lag-0 covariance operator
#' n_grid <- 101
#' 
#' # numerically evaluate lag-h autocovariance operator - WARNING: this takes a minute or so
#' lag <- 0
#' covlagh <- filter_covlagh_operator(theta, lag, n_grid, sigma=sigma)
#' 
#' # # Alternatively simulate with the known eigendecomposition of sigma
#' # covlagh <- filter_covlagh_operator(theta, lag, n_grid, sigma_eigenfunctions = sigma_eigenfunctions, sigma_eigenvalues=sigma_eigenvalues)
#' 
#' # visualise as a surface plot
#' persp(covlagh)
#' 
#' @export
filter_covlagh_operator <- function(theta, lag, n_grid, sigma=NULL, sigma_eigenvalues=NULL, sigma_eigenfunctions=NULL, n_grid_freq=1000){
  
  # if sigma is not defined as kernel, get it from the eigenvalues
  if (is.null(sigma)){
    sigma <- kernel_from_eig( eigenvalues=sigma_eigenvalues, eigenfunctions=sigma_eigenfunctions)
  }
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  int <- array(0, dim=c(n_grid,n_grid))
  sigma_eval <- sigma(grid_matrix,t(grid_matrix))
  for (k in 1:n_grid_freq){
    omega <- pi*k/n_grid_freq # 0..pi
    
    # construct theta discretization as a kernel
    theta_eval <- matrix(0,nrow=n_grid,ncol=n_grid)
    for (ii in 1:n_grid){
      theta_eval[,ii] <- theta(omega, replace(numeric(n_grid), ii, 1))
    }
    
    spec_density <- 1/(2*pi) * theta_eval %*% sigma_eval %*% Conj(t(theta_eval))
    int <- int + spec_density * exp(1i*omega*lag) * pi /n_grid_freq 
    int <- int + t(spec_density) * exp(-1i*omega*lag) * pi /n_grid_freq # contribution on (pi,2pi)
  }
  
  return( Re(int) )
  
}