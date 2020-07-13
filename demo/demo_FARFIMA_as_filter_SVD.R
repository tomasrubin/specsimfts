########################################################################################
# This example demonstrates how to simulate FARFIMA(1,d,0) process whose
# AR operator admits special structure for quick inversion
########################################################################################

library(specsimfts)

###############################################################################################
## set up the FTS dynamics


# innovation covariance operator
sigma <- function(x,y) { pmin(x,y)}

# define filter (frequency response function)
fractional_d <- 0.2
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
fts_x <- filter_simulate(theta, t_max, n_grid, seed_number=seed, sigma=sigma)
end_time <- Sys.time()
print(end_time - start_time)

# plot one trajectory
plot(fts_x[,1], type='l')

## compare the empirical and theoretical autocovariance operator
lag <- 1 # user input here

# calculate the empirical covariance operator
covlagh_empiric <- cov( t(fts_x[,(1+lag):t_max]), t(fts_x[,1:(t_max-lag)]))
persp(covlagh_empiric, ticktype = "detailed") # surface plot. Warning: the visualisation takes long if "n_grid" is high

# theoretical empirical covariance:
covlag0 <- filter_covlagh_operator(theta, 0, n_grid, sigma=sigma)
if (lag == 0){
  covlagh <- covlag0
} else {
  covlagh <- filter_covlagh_operator(theta, lag, n_grid, sigma=sigma)
}

# surface plot of the theoretical covariance
persp(covlagh, ticktype = "detailed") #  Warning: the visualisation takes long if "n_grid" is high

# calculate relative simulation error (in nuclear norm) for this one sample
r <- covlagh_empiric - covlagh
print(paste("Nuclear norm relative error:", sum(svd(r, nu=0, nv=0)$d) / sum(diag(covlag0))))




