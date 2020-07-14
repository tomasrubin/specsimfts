#' Numerically calculate the lag-h covariance operators for functional time series dynamics defined directly by its spectral density operator.
#' The calculation is done by numerically integrating the inverse formula, i.e. the spectral density multiplied by exp(-1i*lag*omega)
#'
#' @title Calculate the lag-h autocovariance operator for FTS defined by its spectral density operator
#' @param spec_density The spectral density operator defined as an integral operator through with the given kernel. Function of three variables, \code{omega}, \code{x}, \code{y}, that assigns the value of the spectral density kernel at point (\code{x},\code{y}) in [0,1]^2 at frequency \code{omega} in [0,pi]. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi].
#' @param lag The lag of the autocovariance to evaluate.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1]^2 for the operator kernel to evaluate.
#' @param n_grid_freq The grid points for the spectral density to evaluate at. Partition of [0,pi].
#' @return lag-h autocovariance operator, matrix of size (\code{n_grid},\code{n_grid})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{CKL_simulate}}
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


spec_density_covlagh_operator <- function(spec_density, lag, n_grid, n_grid_freq=2000){
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  int <- array(0, dim=c(n_grid,n_grid))
  for (k in 1:n_grid_freq){
    omega <- pi*k/n_grid_freq # 0..pi
    spec_density_eval <- spec_density(omega,grid_matrix,t(grid_matrix))
    int <- int + spec_density_eval * exp(1i*omega*lag) *pi /n_grid_freq 
    int <- int + t(spec_density_eval) * exp(-1i*omega*lag) *pi /n_grid_freq 
  }
  
  return( Re(int) )
  
}