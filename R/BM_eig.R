#' This function is not meant to be used by the end-user. It internally provides the decomposition of the Brownian motion for \code{\link{FARFIMA_simulate}} to be run with the Brownian motion covariance.
#'
#' @title Brownian motion kernel eigendecomposition
#' @param n_grid The spatial resolution of the discretization of [0,1], equidistant grid.
#' @return Returns a list of two elements, \code{vectors} is a matrix of size (\code{n_grid},\code{n_grid}) containing the eigenfunctions in the columns of this matrix, and \code{values} vector of eigenvalues.
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{FARFIMA_simulate}}
#' @export
BM_eig <- function(n_grid){
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  
  sigma_eig <- list()
  sigma_eig$vectors <- matrix(0,nrow=n_grid,ncol=n_grid)
  sigma_eig$values <- numeric(n_grid)
  for (n in 1:n_grid){
    sigma_eig$vectors[,n] <- sqrt(2) * sin( (n-0.5)*pi*grid )
    sigma_eig$values[n]  <- 1/(((n-0.5)*pi)^2) 
  }
  
  return(sigma_eig)
}
