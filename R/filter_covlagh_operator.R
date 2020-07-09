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
#' # Define the spectral density operator as an integral operator with kernel
#' k_bbridge <- function(x,y) { pmin(x,y)-x*y }
#' spec_density <- function( omega, x,y ){ 1/(1-0.9 *cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }
#' 
#' # evaluation setting
#' lag <- 1 # change here to evaluate different lag-h autocovariance operator. put "lag <- 0" for lag-0 covariance operator
#' n_grid <- 101
#' 
#' # calculate the lag-h autocovariance operator
#' covlagh <- spec_density_covlagh_operator(spec_density, lag, n_grid)
#' 
#' # visualise as a surface plot
#' persp(covlagh)
#' 
#' @export
filter_covlagh_operator <- function(theta, lag, n_grid, sigma=NULL, sigma_eigenvalues=NULL, sigma_eigenfunctions=NULL, n_grid_freq=500){
  
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
    int <- int + spec_density * exp(1i*omega*lag) *pi /n_grid_freq 
  }
  
  return( Re(int + t(int)) )
  
}