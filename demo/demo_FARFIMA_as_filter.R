###############################################################################################
# This example demonstrates how to simulate FARFIMA(1,d,0) process whose
# AR operator admits special structure for quick inversion
###############################################################################################

library(specsimfts)
## set up the FTS dynamics

# innovation covariance operator (Brownian motion)
sigma_eigenvalues <- function(n) { 1/((n-0.5)*pi)^2 }
sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin((n-0.5)*pi*x) }

# # innovation covariance operator (Brownian bridge)
# sigma_eigenvalues <- function(n) { 1/(n*pi)^2 }
# sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin(n*pi*x) }

# innovation covariance, low rank specification (as list)
# sigma_eigenvalues <- c(1, 0.6, 0.3, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)
# sigma_eigenfunctions <- list(
#   function(x){ sin(2*pi*x) },
#   function(x){ cos(2*pi*x) },
#   function(x){ sin(4*pi*x) },
#   function(x){ cos(4*pi*x) },
#   function(x){ sin(6*pi*x) },
#   function(x){ cos(6*pi*x) },
#   function(x){ sin(8*pi*x) },
#   function(x){ cos(8*pi*x) },
#   function(x){ sin(10*pi*x) },
#   function(x){ cos(10*pi*x) }
# )


# define filter
fractional_d <- +0.2
theta <- function(omega,f){
  ( 2 * sin(omega/2) )^(-fractional_d) *
    (f + (exp(-1i*omega)*0.34) /(1-exp(-1i*omega)*0.34*sqrt(pi)/2*pracma::erfi(1)) *
       rank_one_tensor( function(x){exp((x^2)/2)}, function(x){exp((x^2)/2)}, f ))
}



## simulation setting
t_max <- 1000
n_grid <- 101

# set up random seed
seed <- NULL # no rng seed is inicialized
#seed <- 123

##################################################################################################################
## simulate



## simulate in the spectral domain
start_time <- Sys.time()
fts_x <- filter_simulate(theta, t_max, n_grid, seed_number=seed, sigma_eigenfunctions = sigma_eigenfunctions, sigma_eigenvalues=sigma_eigenvalues)
end_time <- Sys.time()
print(end_time - start_time)

# display the first curve
plot(fts_x[,1], type='l')

## compare the empirical and theoretical autocovariance operator
lag <- 0 # user input here

# calculate the empirical covariance operator
covlagh_empiric <- cov( t(fts_x[,(1+lag):t_max]), t(fts_x[,1:(t_max-lag)]))
persp(covlagh_empiric, ticktype = "detailed") # surface plot. Warning: the visualisation takes long if "n_grid" is high

# theoretical empirical covariance - WARNING: this cat take a minute or longer !!!
covlag0 <- filter_covlagh_operator(theta, 0, n_grid, sigma_eigenvalues=sigma_eigenvalues, sigma_eigenfunctions=sigma_eigenfunctions)
if (lag == 0){
  covlagh <- covlag0
} else {
  covlagh <- filter_covlagh_operator(theta, lag, n_grid, sigma_eigenvalues=sigma_eigenvalues, sigma_eigenfunctions=sigma_eigenfunctions)
}

# surface plot of the theoretical covariance
persp(covlagh, ticktype = "detailed") #  Warning: the visualisation takes long if "n_grid" is high

# calculate relative simulation error (in nuclear norm) for this one sample
r <- covlagh_empiric - covlagh
print(paste("Nuclear norm relative error:", sum(svd(r, nu=0, nv=0)$d) / sum(diag(covlag0))))




