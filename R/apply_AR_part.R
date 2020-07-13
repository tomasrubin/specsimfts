#' This function is not meant to be used by the end-user. It is internal function for the runtime of the hybrid simulation withing the FARFIMA_simulate function.
#' 
#'
#' @title Apply the autoregressive recursion in the temporal domain
#' @param fts_x The functional time series, represented as a matrix of size (\code{n_grid},\code{t_max}), which the autoregressive recursion is applied to. 
#' @param operators_ar The list of length 'p' the order of the autoregressive part. The autoregressive operators are considered to be integral operators defined through their kernels which are saved as the elements of the list \code{operators_ar} as functions of two variables, \code{x} and \code{y}, returning the value of the kernel at point (\code{x},\code{y}). In case of degenerate autoregressive part define \code{operators_ar} as an empty list.
#' @return functional time series sample, matrix of size (\code{n_grid},\code{t_max}) where these sample size parameters are copied from the input parameter \code{fts_x}
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{FARFIMA_simulate}}
#' @export
apply_AR_part <- function(fts_x, operators_ar){
  
  # first "ar_order" elements of the ouput FTS are zeros
  
  # grid for evaluation
  n_grid <- nrow(fts_x)
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # shortcut for AR oder
  ar_order <- length(operators_ar)
  
  # do the AR recursion only if nondegenerate AR part
  if (ar_order > 0){
    
    t_max <- ncol(fts_x)
    fts_x_new <- matrix(0, nrow=nrow(fts_x), ncol=t_max)
    
    # express AR operators
    ar_evaluation <- operators_ar
    for (j in 1:ar_order){
      ar_evaluation[[j]] <- operators_ar[[j]](grid_matrix,t(grid_matrix)) 
    }
    
    # AR cycle
    for (t in (1+ar_order):t_max){
      fts_x_new[,t] <- fts_x[,t] # noise contribution
      for (j in 1:ar_order){
        fts_x_new[,t] <- fts_x_new[,t] + ar_evaluation[[j]] %*% fts_x_new[,t-j] / n_grid
      }
    }
    
    return(fts_x_new)
    
  } else {
    return(fts_x)
  }
  
}