###############################################################################################
# This example demonstrates how to simulate a filtered white noise process defined
# using a custom filter (definition in the comments below) and Brownian motion innovation covariance.
###############################################################################################

library(specsimfts)
## set up the FTS dynamics

# innovation covariance operator (Brownian motion)
sigma <- function(x,y) { pmin(x,y) }

# # innovation covariance operator (Brownian bridge)
# sigma <- function(x,y) { pmin(x,y) - x*y }


# define filter by its frequency response function
# (Theta(\omega) f)(x) = 2f(x) + i*f(1-x) + omega \int_0^x f(y)dy + ((v_1) \otimes (v_2))(f)(x) + \int_0^1 K_omega(x,y)f(y)dy
# with (v1)(x) = sin(x), (v2)(x) = exp(x), andK_omega(x,y) = sin(omega+x+2y)
theta <- function(omega,f){
  2*f+ # 2f(x)
    1i*rev(f) + # i*f(1-x)
    omega*cumsum(f)/length(f) + # omega * \int_0^x f(y)dy
    rank_one_tensor( function(x){sin(x)}, function(x){exp(x)}, f ) + # ((v_1) \otimes (v_2))(f)(x) with (v1)(x) = sin(x) and (v2)(x) = exp(x)
    kernel_operator( function(x,y){sin(omega+x+2*y)}, f ) # \int_0^1 K(x,y)f(y)dy with K_omega(x,y) = sin(omega+x+2y)
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

# display the first curve
plot(fts_x[,1], type='l')

## compare the empirical and theoretical autocovariance operator
lag <- 0 # user input here

# calculate the empirical covariance operator - the numerical integration can take long !!!
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




