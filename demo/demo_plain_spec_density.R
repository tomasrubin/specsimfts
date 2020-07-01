# This demo file shows simulation of the functional time series who's
# spectral density kernels given by a direct formula


## Define the functions for harmonic eigenvalues and eigenfunctions
k_bbridge <- function(x,y) { pmin(x,y)-x*y }
spec_density <- function( omega, x,y ){ 1/(1-0.9 *cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }

## simulation setting
t_max <- 1000 # time horizon to be simulated
n_grid <- 101 # spatial resolution for visualisation on discretisation of [0,1]. warning: scales badly with high "n_grid"
n_pc <- n_grid # number of numerically calculated eigenvalues to use. there is negligible computational gain, thus "n_pc = n_grid" is recommended

# set up random seed
seed <- NULL # no random seed inicialized
#seed <- 123


########################################################################################
## simulate trajectory
start_time <- Sys.time()
fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc, seed_number = seed)
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
covlag0 <- spec_density_covlagh_operator(spec_density, 0, n_grid)
if (lag == 0){
  covlagh <- covlag0
} else {
  covlagh <- spec_density_covlagh_operator(spec_density, lag, n_grid)
}

# surface plot of the theoretical covariance
persp(covlagh, ticktype = "detailed") #  Warning: the visualisation takes long if "n_grid" is high

# calculate relative simulation error (in nuclear norm) for this one sample
r <- covlagh_empiric - covlagh
print(paste("Nuclear norm relative error:", sum(svd(r, nu=0, nv=0)$d) / sum(diag(covlag0))))

