#' Numerically calculate the lag-h covariance operators for functional time series dynamics defined by its harmonic Karhunen-Loeve expansion.
#' The calculation is done by numerically integrating the inverse formula, i.e. the spectral density multiplied by exp(-1i*lag*omega)
#'
#' @title Calculate the lag-h autocovariance operator for functional time series dynamics defined by its harmonic Karhunen-Loeve expansion.
#' @param harmonic_eigenvalues function of two variables, \code{omega} and \code{n}, that assigns the \code{n}-th harmonic eigenvalue at frequency \code{omega}. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi]..
#' @param harmonic_eigenfunctions function of three variables, \code{omega}, \code{n} and \code{x}, that assigns the \code{n}-th harmonic eigenfunction at point \code{x} in [0,1] at frequency \code{omega}. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi].
#' @param lag The lag of the autocovariance to evaluate.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1]^2 for the operator kernel to evaluate.
#' @param n_pc The number of harmonic eigenfunctions to be used for the numerical integration at each frequency.
#' @param n_grid_freq The grid points for the spectral density to evaluate at. Partition of [0,pi].
#' @return lag-h autocovariance operator, matrix of size (\code{n_grid},\code{n_grid})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{HKL_simulate}}
#' @examples
#' # Define the eigenvalues and eigenfunction of the functional time series
#' harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
#' harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*x-omega)  ) }
#' 
#' # evaluation setting
#' lag <- 1 # change here to evaluate different lag-h autocovariance operator. put "lag <- 0" for lag-0 covariance operator
#' n_grid <- 101
#' n_pc <- 100
#' 
#' # calculate the lag-h autocovariance operator
#' covlagh <- HKL_covlagh_operator(harmonic_eigenvalues, harmonic_eigenfunctions, lag, n_grid, n_pc)
#' 
#' # visualise as a surface plot
#' persp(covlagh)
#' 
#' @export


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