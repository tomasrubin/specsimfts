% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FARFIMA_covlagh_operator.R
\name{FARFIMA_covlagh_operator}
\alias{FARFIMA_covlagh_operator}
\title{Calculate the lag-h autocovariance operator of the FARFIMA(p,d,q) process}
\usage{
FARFIMA_covlagh_operator(FARFIMA_pars, lag, n_grid, n_grid_freq = 500)
}
\arguments{
\item{FARFIMA_pars}{The list of the parameters for the FARFIMA(p,d,q) process. Must contain fields: (i) \code{fractional_d}, a real number in the open interval (-0.5,0.5) controling the fractional integration degree. \code{fractional_d} being positive corresponds to long-rande dependence behaviour. (ii) \code{operators_ar}, the list of length 'p' the order of the autoregressive part. The autoregressive operators are considered to be integral operators defined through their kernels which are saved as the elements of the list \code{operators_ar} as functions of two variables, \code{x} and \code{y}, returning the value of the kernel at point (\code{x},\code{y}). In case of degenerate autoregressive part define \code{operators_ar} as an empty list. (iii) \code{operators_ma}, the list of length 'q', the order of the moving average part. Just like \code{operators_ar} its a liks of functions - the kernels of the moving average operators. (iv) The covariance opperator of the stochastic innovation process can be defined either through (iv-a) its kernel,  (iv-b) finite rank eigendecomposition, (iv-c) infinite rank decomposition. In the case (iv-a), define \code{sigma} as a function of two variables \code{x} and \code{y}, returning the value of the covariance kernel at point (\code{x},\code{y}). In the case (iv-b), define the elements \code{sigma_eigenvalues} as a vector of finitely many eigenvalues and \code{sigma_eigenfunctions} as a list of the same length as \code{sigma_eigenvalues} with each element being a function of variable \code{x} returning the value of that eigenfunction at point \code{x}. In the case (iv-c), define the elements \code{sigma_eigenvalues} as a function of the variable \code{n} returning the \code{n}-th eigenvalue and the element \code{sigma_eigenfunctions} as a function of two variables, \code{n} and \code{x}, returning the value of the \code{n}-th eigenfunctions at point \code{x}. See the example bellow for some examples on how to set up \code{FARFIMA_pars}.}

\item{lag}{The lag of the autocovariance to evaluate.}

\item{n_grid}{Number of grid points (spatial resolution) of the discretisation of [0,1]^2 for the operator kernel to evaluate.}

\item{n_grid_freq}{The grid points for the spectral density to evaluate at. Partition of [0,pi].}
}
\value{
lag-h autocovariance operator, matrix of size (n_grid,n_grid)
}
\description{
Numerically calculate the lag-h covariance operators for FARFIMA(p,d,q) process. The calculation is done by numerically integrating the inverse formula, i.e. the spectral density multiplied by exp(-1i*lag*omega). If the process has non-degenerate autoregressive part, the evaluation of the spectral density requires matrix inversion at each frequency. The function is very slow for large or even moderate \code{n_grid}.
}
\examples{
# (i) fractional integration
fractional_d <- 0 # in the open interval (-0.5, 0.5), positive number means long-range dependence

# (ii) autoregressive operators
operators_ar <- list(
function(x,y){ 0.3*sin(x-y) },
function(x,y){ 0.3*cos(x-y) },
function(x,y){ 0.3*sin(2*x) },
function(x,y){ 0.3*cos(y) }
)
# operators_ar <- list() # use empty list for degenerate AR part

# (iii) moving average kernels
# you can put here arbitrary long list of operators
operators_ma <- list(
function(x,y){ x+y },
function(x,y){ x },
function(x,y){ y }
)
# operators_ma <- list() # use empty list for degenerate MA part

# (iv-b) covariance of the inovation defined through eigenvalues and eigenfunctions
# you can put here arbitrary long lists but their lenghts should match
sigma_eigenvalues <- c(1, 0.6, 0.3, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)
sigma_eigenfunctions <- list(
function(x){ sin(2*pi*x) },
function(x){ cos(2*pi*x) },
function(x){ sin(4*pi*x) },
function(x){ cos(4*pi*x) },
function(x){ sin(6*pi*x) },
function(x){ cos(6*pi*x) },
function(x){ sin(8*pi*x) },
function(x){ cos(8*pi*x) },
function(x){ sin(10*pi*x) },
function(x){ cos(10*pi*x) }
)

# # (iv-c) innovation covariance operator (Brownian motion)
# sigma_eigenvalues <- function(n) { 1/((n-0.5)*pi)^2 }
# sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin((n-0.5)*pi*x) }

# put the parameters into one list
FARFIMA_pars <- list(fractional_d=fractional_d, operators_ar=operators_ar, operators_ma=operators_ma, sigma_eigenvalues=sigma_eigenvalues,sigma_eigenfunctions=sigma_eigenfunctions)

# # (iv-a) Alternatively, define the kernel of the white noise innovation.
# sigma <- function(x,y) { pmin(x,y) } # Brownian motion
# FARFIMA_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma, sigma=sigma)


if (FARFIMA_test_stationarity(FARFIMA_pars)){
 # calculate the lag-h autocovariance kernel
 lag <- 1 # change here to evaluate different lag-h autocovariance operator. put "lag <- 0" for lag-0 covariance operator
 covlagh <- FARFIMA_covlagh_operator(FARFIMA_pars, lag, n_grid)
   
 # surface plot the covlagh
 persp(covlagh)
}
 

}
\references{
Rubin, Panaretos. \emph{Simulation of stationary functional time series with given spectral density}. arXiv, 2020
}
\seealso{
\code{\link{FARFIMA_simulate}}, \code{\link{FARFIMA_test_stationarity}}
}
