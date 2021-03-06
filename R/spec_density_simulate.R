#' Simulate functional time series sample defined through their spectral density. The simulation routine needs to
#' perform the SVD decomposition at each freaquency, therefore the simulation is quite slow. For a speedy simulation
#' try to define the spectral density operator through its harmonic Karhunen-Loeve expansion (\code{\link{spec_density_simulate}}), or by the
#' white noise filter approach (\code{\link{spec_density_simulate}}).
#' 
#'
#' @title Simulate a FTS given by its spectral density operator
#' @param spec_density The spectral density operator defined as an integral operator through with the given kernel. Function of three variables, \code{omega}, \code{x}, \code{y}, that assigns the value of the spectral density kernel at point (\code{x},\code{y}) in [0,1]^2 at frequency \code{omega} in [0,pi]. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi].
#' @param t_max Time horizon to be simulated. Must be an even number, otherwise it is increased by one.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1] where the FTS is to be simulated.
#' @param n_pc The number of harmonic eigenfunctions to be used for the simulation at each frequency.
#' @param seed_number The random seed inicialization for the simulation. The value \code{NULL} means no inicialization
#' @param include_freq_zero If set \code{TRUE}, the zero frequency is included for simulation in the spectral domain. Set \code{FALSE} for processes with singularity at frequency zero, e.g. the long-range dependent FARFIMA(p,d,q) process.
#' @return functional time series sample, matrix of size (\code{n_grid},\code{t_max})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{spec_density_covlagh_operator}}, \code{\link{CKL_simulate}}
#' @examples
#' # Define the spectral density operator as an integral operator with kernel
#' k_bbridge <- function(x,y) { pmin(x,y)-x*y }
#' spec_density <- function( omega, x,y ){ 1/(1-0.9 *cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }
#' 
#' # simulation setting
#' t_max <- 1000 # time horizon to be simulated
#' n_grid <- 101 # spatial resolution for visualisation on discretisation of [0,1]. warning: scales badly with high "n_grid"
#' n_pc <- n_grid # number of numerically calculated eigenvalues to use. there is negligible computational gain, thus "n_pc = n_grid" is recommended
#' 
#' # simulate a sample 
#' fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc)
#' 
#' # display the first curve
#' plot( fts_x[,1], type='l' )
#' 
#' @export

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
  for (ii in 1:(t_half+1)){
    
    if (ii == t_half + 1){
      # zero frequency
    } else {
      # normal frequency
      omega <- (2*ii*pi)/t_max
    }
    
    # SVD decomposition of spec density
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