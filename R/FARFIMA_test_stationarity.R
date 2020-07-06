#' Numerically verify if the FARFIMA(p,d,q) process is stationary. The stationarity depends solely on the autoregressive part. The method construct the order-1 autoregressive process in the state space and calculate the 1000-th power of the state-space autoregressive operator. If its norm is less than one, the process is pronnounced stationary.
#'
#' @title Check stationarity of the FARFIMA(p,d,q) process
#' @param FARFIMA_pars The list of the parameters for the FARFIMA(p,d,q) process. Must contain fields: (i) \code{fractional_d}, a real number in the open interval (-0.5,0.5) controling the fractional integration degree. \code{fractional_d} being positive corresponds to long-rande dependence behaviour. (ii) \code{operators_ar}, the list of length 'p' the order of the autoregressive part. The autoregressive operators are considered to be integral operators defined through their kernels which are saved as the elements of the list \code{operators_ar} as functions of two variables, \code{x} and \code{y}, returning the value of the kernel at point (\code{x},\code{y}). In case of degenerate autoregressive part define \code{operators_ar} as an empty list. (iii) \code{operators_ma}, the list of length 'q', the order of the moving average part. Just like \code{operators_ar} its a liks of functions - the kernels of the moving average operators. (iv) The covariance opperator of the stochastic innovation process can be defined either through (iv-a) its kernel,  (iv-b) finite rank eigendecomposition, (iv-c) infinite rank decomposition. In the case (iv-a), define \code{sigma} as a function of two variables \code{x} and \code{y}, returning the value of the covariance kernel at point (\code{x},\code{y}). In the case (iv-b), define the elements \code{sigma_eigenvalues} as a vector of finitely many eigenvalues and \code{sigma_eigenfunctions} as a list of the same length as \code{sigma_eigenvalues} with each element being a function of variable \code{x} returning the value of that eigenfunction at point \code{x}. In the case (iv-c), define the elements \code{sigma_eigenvalues} as a function of the variable \code{n} returning the \code{n}-th eigenvalue and the element \code{sigma_eigenfunctions} as a function of two variables, \code{n} and \code{x}, returning the value of the \code{n}-th eigenfunctions at point \code{x}. See the example bellow for some examples on how to set up \code{FARFIMA_pars}.
#' @param n_grid Number of grid points (spatial resolution) of the discretisation of [0,1]^2. The method checks if the process is stationary on the discretization level. Assuming smoothness, it shoudn't be dependent on the grid resolution unless very coarse.
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @return Returns \code{TRUE} or \code{FALSE} if stationary or not.
#' @seealso \code{\link{FARFIMA_simulate}} \code{\link{FARFIMA_test_stationarity}}
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
#' # test stationarity
#' print(FARFIMA_test_stationarity(FARFIMA_pars))
#'  
#' 
#' @export
FARFIMA_test_stationarity <- function(FARFIMA_pars, n_grid=101){
  
  # evaluate operators on grid
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # save model order
  ar_order <- length(FARFIMA_pars$operators_ar)
  
  
  if (ar_order == 0){# if no AR part, always stationary
    return(TRUE)
  } else {
    
    # express AR operators
    operators_ar_eval <- FARFIMA_pars$operators_ar
    
    for (j in 1:ar_order){
      operators_ar_eval[[j]] <- FARFIMA_pars$operators_ar[[j]](grid_matrix,t(grid_matrix)) / n_grid
    }
    
    # compose the AR operator in the state-space
    composed_a <- matrix(0, nrow=ar_order*n_grid, ncol=ar_order*n_grid)
    for (ii in 1:ar_order){
      composed_a[ 1:n_grid , ((ii-1)*n_grid+1):(ii*n_grid)] <- operators_ar_eval[[ii]]
    }
    if (ar_order>1){
      for (ii in 1:(ar_order-1)){
        composed_a[ (ii*n_grid+1):((ii+1)*n_grid), ((ii-1)*n_grid+1):(ii*n_grid) ] <- diag(n_grid)
      }
    }
    
    m <- composed_a %^% 1000 # perform matrix power
    
    if (is.na(max(m))){
      return(F)
    } else {
      if (is.finite(m[1,1])){
        return( norm( m, type="2" ) < 1 )
      } else {
        return(F)
      }
    }
  }
}