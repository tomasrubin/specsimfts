#' Simulate functional time series sample given its harmonic Karhunen-Loeve decomposition.
#'
#' @title Simulate FTS given its harmonic KL decomposition
#' @param harmonic_eigenvalues function of two variables, \code{omega} and \code{n}, that assigns the \code{n}-th harmonic eigenvalue at frequency \code{omega}. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi]..
#' @param harmonic_eigenfunctions function of three variables, \code{omega}, \code{n} and \code{x}, that assigns the \code{n}-th harmonic eigenfunction at point \code{x} in [0,1] at frequency \code{omega}. Must be well defined for frequencies (0,pi]. The interval [pi,2pi) is not used and is calculated by mirroring of (0,pi].
#' @param t_max Time horizon to be simulated. Must be an even number, otherwise it is increased by one.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1] where the FTS is to be simulated.
#' @param n_pc The number of harmonic eigenfunctions to be used for the simulation at each frequency.
#' @param seed_number The random seed inicialization for the simulation. The value "NULL" means no inicialization
#' @param include_freq_zero If set \code{TRUE}, the zero frequency is included for simulation in the spectral domain. Set \code{FALSE} for processes with singularity at frequency zero, such as the long-range dependent FARFIMA(p,d,q) process.
#' @return functional time series sample, matrix of size (\code{n_grid},\code{t_max})
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{HKL_covlagh_operator}} \code{\link{spec_density_simulate}}
#' @examples
#' # Define the eigenvalues and eigenfunction of the process to simulate
#' harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
#' harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*x-omega)  ) }
#' 
#' # simulation setting
#' t_max <- 1000
#' n_grid <- 101
#' n_pc <- 100
#' 
#' # simulate trajectory
#' fts_x <- f_CKL_simulate(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc)
#' 
#' # display the first curve
#' plot( fts_x[,1], type='l' )
#' 
#' @export
HKL_simulate <- function(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc, seed_number = NULL, include_freq_zero = F){
  
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
  
  # get Vs
  vs <- matrix(0, n_grid, t_max)
  #pb <- txtProgressBar(style = 3)
  for (n in (1:n_pc)){
    #setTxtProgressBar(pb, n/n_pc)
    # prepare eigenvalue and eigenfunction of Sigma
    
    for (ii in 1:(t_half+1)){
      
      # split by if we're at zero freq or not
      if (ii == t_half + 1){
        # this ii corresponds to t_max = zero frequency
        if (include_freq_zero){
          omega <- 0
          phi <- harmonic_eigenfunctions(omega, n, grid )
          vs[,t_max] <- vs[,t_max] + sqrt(harmonic_eigenvalues(omega,n)) * phi * zs[n,ii]
        }
        # do nothing
      } else {
        # all the other ii, going from 1 to t_half
        
        # function
        omega <- (2*ii*pi)/t_max
        phi <- harmonic_eigenfunctions(omega, n, grid )
        
        vs[,ii] <- vs[,ii] + sqrt(harmonic_eigenvalues(omega,n)) * phi * zs[n,ii]
        
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
