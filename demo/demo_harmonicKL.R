# This demo file shows simulation of the functional time series who's
# spectral density operators are given by their eigendecomposition.


## Define the functions for harmonic eigenvalues and eigenfunctions
harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*x-omega)  ) }

## simulation setting
t_max <- 1000 # time horizon to be simulated
n_grid <- 101 # spatial resolution for visualisation on discretisation of [0,1]
n_pc <- 100 # number of harmonic principal components (eigenfunctions) to use for simulation

########################################################################################
## simulate trajectory
start_time <- Sys.time()
fts_x <- HKL_simulate(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc)
end_time <- Sys.time()
print( difftime(end_time,start_time, units="secs")) # print the time difference

## display the first curve
par(mfrow=c(1,3))
plot( fts_x[,1], type='l' )

#########################################################################################
## compare the empirical and theoretical autocovariance operator
lag <- 1 # user input here

# calculate the empirical covariance operator
covlagh_empiric <- cov( t(fts_x[,(1+lag):t_max]), t(fts_x[,1:(t_max-lag)]))
persp(covlagh_empiric, ticktype = "detailed") # surface plot. Warning: the visualisation takes long if "n_grid" is high

# theoretical empirical covariance
covlag0 <- HKL_covlagh_operator(harmonic_eigenvalues, harmonic_eigenfunctions, 0, n_grid, n_pc=100) # you may edit n_pc for more precice evaluation
if (lag == 0){
  covlagh <- covlag0
} else {
  covlagh <- HKL_covlagh_operator(harmonic_eigenvalues, harmonic_eigenfunctions, lag, n_grid, n_pc=100) # you may edit n_pc for more precice evaluation
}

# surface plot of the theoretical covariance
persp(covlagh, ticktype = "detailed") #  Warning: the visualisation takes long if "n_grid" is high

# calculate relative simulation error (in nuclear norm) for this one sample
r <- covlagh_empiric - covlagh
print(paste("Nuclear norm relative error:", sum(svd(r, nu=0, nv=0)$d) / sum(diag(covlag0))))

