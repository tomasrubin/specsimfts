#' This function applies the integral operator with the kernel \code{ker} onto the function \code{f}. Can be useful for defining custom filter to be used for the filtered white noise simulation approach, see \code{\link{filter_simulate}}.
#' 
#'
#' @title Apply the integral operator defined by its kernel
#' @param ker The kernel function of the integral operator defined as a function of two variables, \code{x} and \code{y}, returning the value of the integral kernel at the point (\code{x},\code{y}). It is evaluated at the grid of [0,1]^2 where the resolution is matched to the length of the vector \code{f}.
#' @param f The discretized function \code{f} represented as a vector. Assumes equidistant partition of [0,1]. 
#' @return Vector of the same size as \code{f} representing the discretized function - the result of the application "ker(f)".
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{rank_one_tensor}}, \code{\link{filter_simulate}}
#' @examples 
#' # define kernel
#' ker <- function(x,y){ exp(x^2+y^2) }
#' 
#' # discretize sinus function
#' grid <- seq(0,1, length.out = 101)
#' f <- sin(grid)
#' 
#' # apply onto function
#' plot(f, type="l", col="blue", ylim=c(-0.5,2.5))
#' lines( kernel_operator(ker,f), col="red" )
#' legend("topleft",c("f","ker(f)"), col=c("blue","red"), lty=1)
#' 
#' # Note that the above kernel function is in fact a rank one tensor, thus the application could be implemented as (which is faster to evaluate): rank_one_tensor( function(x) exp(x^2), function(x) exp(x^2), f )
#' @export
kernel_operator <- function( ker, f ){
  grid <- seq(0,1,length.out=length(f))
  grid_matrix <- kronecker(grid,matrix(1,1,length(f)))
  ker( grid_matrix, t(grid_matrix) ) %*% f / length(grid)
}