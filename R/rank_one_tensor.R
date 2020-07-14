#' This function applies the operator given as a rank-one tensor defines rank one operator defined as: "(g1 'otimes' g2)(f) = <f,g2>g1"
#' 
#'
#' @title Apply the rank-one-tensor onto a function
#' @param g1 The function \code{g1} for the rank-one-tensor define as a function of variable \code{x} returning the the value "g1(x)".
#' @param g2 The function \code{g2} for the rank-one-tensor define as a function of variable \code{x} returning the the value "g2(x)".
#' @param f The discretized function \code{f} represented as a vector. Assumes equidistant partition of [0,1]. 
#' @return Vector of the same size as \code{f} representing the discretized function - the result of the application "(g1 'otimes' g2)(f)".
#' @references Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
#' @seealso \code{\link{kernel_operator}} \code{\link{filter_simulate}}
#' @examples 
#' # define kernel
#' g1 <- function(x){ exp(x^2) }
#' g2 <- function(x){ sin(x) }
#' 
#' # discretize sinus function
#' grid <- seq(0,1, length.out = 101)
#' f <- sin(grid)
#' 
#' # apply onto function
#' plot(f, type="l", col="blue", ylim=c(-0.5,1))
#' lines( rank_one_tensor(g1,g2,f), col="red" )
#' legend("topleft",c("f","<f,g2>g1"), col=c("blue","red"), lty=1)
#' @export
rank_one_tensor <- function(g1, g2, f){
  grid <- seq(0,1,length.out=length(f))
  return( mean( f * g2(grid) ) * g1(grid) )
}
