# This example demonstrates how to simulate FARFIMA(p,d,q) process whose
# MA and AR parameters are defined as integral operators through their kernels
# and innovation noise covariance is given by its eigendecomposition

###############################################################################################
## set up the FTS dynamics

# fractional integration
fractional_d <- 0 # in the open interval (-0.5, 0.5), positive number means long-range dependence

# autoregressive operators
operators_ar <- list(
  function(x,y){ 0.3*sin(x-y) },
  function(x,y){ 0.3*cos(x-y) },
  function(x,y){ 0.3*sin(2*x) },
  function(x,y){ 0.3*cos(y) }
)
# operators_ar <- list() # use empty list for degenerate AR part

# moving average kernels
# you can put here arbitrary long list of operators
operators_ma <- list(
  function(x,y){ x+y },
  function(x,y){ x },
  function(x,y){ y }
)
# operators_ma <- list() # use empty list for degenerate MA part

# covariance of the inovation defined through eigenvalues and eigenfunctions
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


## alternatively, you can define the eigenvalues as functions of "n" and eigenfunctions as bivariate functions of "n" and "x"
## uncomment one of the two following kernels and comment the low-rank eigendecomposition above

# # innovation covariance operator (Brownian motion)
# sigma_eigenvalues <- function(n) { 1/((n-0.5)*pi)^2 }
# sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin((n-0.5)*pi*x) }

# # innovation covariance operator (Brownian bridge)
# sigma_eigenvalues <- function(n) { 1/(n*pi)^2 }
# sigma_eigenfunctions <- function(n,x) { sqrt(2)*sin(n*pi*x) }


# put the parameters into one list
FARFIMA_pars <- list(fractional_d=fractional_d,
                     operators_ar=operators_ar,
                     operators_ma=operators_ma,
                     sigma_eigenvalues=sigma_eigenvalues,
                     sigma_eigenfunctions=sigma_eigenfunctions)


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
  start_time <- Sys.time()
  fts_x <- FARFIMA_simulate(FARFIMA_pars, t_max, n_grid, seed_number=seed, hybrid_ar = F)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # plot one trajectory
  par(mfrow=c(1,3))
  plot(fts_x[,1], type='l')
  
  ## compare the empirical and theoretical autocovariance operator
  lag <- 1 # user input here
  
  # calculate the empirical covariance operator
  covlagh_empiric <- cov( t(fts_x[,(1+lag):t_max]), t(fts_x[,1:(t_max-lag)]))
  persp(covlagh_empiric, ticktype = "detailed") # surface plot. Warning: the visualisation takes long if "n_grid" is high
  
  # theoretical empirical covariance:
  # WARNING: the evaluation of true autocovariance operators is !!extremely slow!! for high "n_grid" !!!
  covlag0 <- FARFIMA_covlagh_operator(FARFIMA_pars, 0, n_grid)
  if (lag == 0){
    covlagh <- covlag0
  } else {
    covlagh <- FARFIMA_covlagh_operator(FARFIMA_pars, lag, n_grid)
  }
  
  # surface plot of the theoretical covariance
  persp(covlagh, ticktype = "detailed") #  Warning: the visualisation takes long if "n_grid" is high
  
  # calculate relative simulation error (in nuclear norm) for this one sample
  r <- covlagh_empiric - covlagh
  print(paste("Nuclear norm relative error:", sum(svd(r, nu=0, nv=0)$d) / sum(diag(covlag0))))

} else {
  print("Not stationary AR part.")
}


