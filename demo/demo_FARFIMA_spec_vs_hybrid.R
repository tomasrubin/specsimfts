########################################################################################
# This example demonstrates how to simulate FARFIMA(p,d,q) process
# using the fully spectral and hybrid regularisation
########################################################################################

library(specsimfts)

###############################################################################################
## set up the FTS dynamics

# fractional integration
fractional_d <- 0.3 # in the open interval (-0.5, 0.5), positive number means long-range dependence

# autoregressive operators
operators_ar <- list(
  function(x,y){ 0.3*exp(x+y) }
)
# operators_ar <- list() # use empty list for degenerate AR part

# moving average kernels
# you can put here arbitrary long list of operators
operators_ma <- list(
  function(x,y){ x+y }
)
# operators_ma <- list() # use empty list for degenerate MA part

# covariance of the inovation defined by the kernel
# sigma <- function(x,y) { exp( -(x-y)^2 ) } # exponential kernel
sigma <- function(x,y) { pmin(x,y) } # Brownian motion
# sigma <- function(x,y) { pmin(x,y) - x*y } # Brownian bridge

# put the parameters into one list
FARFIMA_pars <- list(fractional_d=fractional_d,
                     operators_ar=operators_ar,
                     operators_ma=operators_ma,
                     sigma=sigma)


## set simulation setting
t_max <- 1600
n_grid <- 101

# set up random seed
seed <- NULL # no rng seed is inicialized
#seed <- 123

##################################################################################################################
## simulate

# first test stationarity of the defined autoregressive part
if (FARFIMA_test_stationarity(FARFIMA_pars, 101)){

  ## simulate in the spectral domain
  spec_start_time <- Sys.time()
  spec_fts_x <- FARFIMA_simulate(FARFIMA_pars, t_max, n_grid, seed_number=seed, hybrid_ar = F)
  spec_end_time <- Sys.time()
  
  ## simulate by the hybrid method
  hybrid_start_time <- Sys.time()
  hybrid_fts_x <- FARFIMA_simulate(FARFIMA_pars, t_max, n_grid, seed_number=seed, hybrid_ar = T)
  hybrid_end_time <- Sys.time()
  
  ## uncomment these to display trajectories
  # as the methods are conceptually different, the trajectories are different even if inicialized by the same seed
  # plot(spec_fts_x[,1], type='l', main="fully spectral")
  # plot(hybrid_fts_x[,1], type='l', main="hybrid")
  
  ## display timing
  print(paste("Fully spectral (solving lin.eq. at each frequency):",round(spec_end_time - spec_start_time,2),"seconds."))
  print(paste("Hybrid simulation:",round(hybrid_end_time - hybrid_start_time,2),"seconds."))
  
} else {
  print("Not stationary AR part.")
}


