###################################################################################
## simulate FARIMA(1,d,0) using either yes or no the analytic decomposition of Brownian motion covariance, and either yes or no rank-one inversion of the AR operator
# how_to_AR = "sm" ...... Sherman-Morrison trick
# how_to_AR = "spec" .... fully spectral
# how_to_AR = "hybrid" .. apply the AR recursion in the temporal domain
simulation_run_FARIMA_PC <- function( analytic_eig, how_to_AR, t_max, n_grid, lag_to_compare, seed_number=NaN ){
  
  # here I'm saving empirical col kernels for the designated lags, for each sim run
  lag_to_compare_n <- length(lag_to_compare)
  covlags_all <- array(NaN, dim=c(lag_to_compare_n,n_grid,n_grid) )
  
  # simulate trajectory
  start_time = Sys.time()
  sigma <- function(x,y) { pmin(x,y)}
  if (how_to_AR == "sm"){
    
    # define filter
    fractional_d <- 0.2
    theta <- function(omega,f){
      ( 2 * sin(omega/2) )^(-fractional_d) *
        (f + (exp(-1i*omega)*0.34) /(1-exp(-1i*omega)*0.34*sqrt(pi)/2*erfi(1)) *
           rank_one_tensor( function(x){exp((x^2)/2)}, function(x){exp((x^2)/2)}, f ))
    }
    
    if (analytic_eig){
      fts_x <- filter_simulate(sigma, theta, t_max, n_grid, seed_number=seed_number, sigma_eig = BM_eig(n_grid), include_zero_freq=F)
    } else {
      fts_x <- filter_simulate(sigma, theta, t_max, n_grid, seed_number=seed_number, include_zero_freq=F)
    }
    
    
  } else {
    
     # how_to_AR == "spec" or "hybrid"
      
      # simulate FARFIMA(0,0.2,0)
      fractional_d <- 0.2
      operators_ar <- list(
        function(x,y){ -0.34*exp( (x^2+y^2)/2 ) }
      )
      operators_ma <- list()
      if (analytic_eig){
        farima_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma,sigma="bm")
      } else {
        farima_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma,sigma=sigma)
      }

      fts_x <- FARIMA_simulate(farima_pars, t_max, n_grid, seed_number=seed, hybrid_ar = (how_to_AR == "hybrid"))
      
  }
  
  end_time  = Sys.time()
  timing <- difftime(end_time,start_time, units="secs")  # saving the simulation time
  
  # evaluate and save empirical cov kernels
  for (lag_i in 1:lag_to_compare_n){
    lag <- lag_to_compare[lag_i]
    
    # calculate covariance kernel
    m <- cov( t(fts_x[,(lag+1):t_max]), t(fts_x[,1:(t_max-lag)]) )
    covlags_all[lag_i,,] <- m
  }
  
  return(list(covlags_all=covlags_all,timing=timing))
}


###################################################################################
## simulate ARMA(4,3)
# method_i = 1 ... lowrank_spec
# method_i = 2 ... lowrank_hybrid
# method_i = 3 ... svd_spec
# method_i = 4 ... svd_hybrid
simulation_run_ARMA_lowrank <- function( t_max, n_grid, lag_to_compare, method_i, seed_number=NaN ){
  
  # here I'm saving empirical col kernels for the designated lags, for each sim run
  lag_to_compare_n <- length(lag_to_compare)
  covlags_all <- array(NaN, dim=c(lag_to_compare_n,n_grid,n_grid) )
  
  # simulate trajectory
  start_time = Sys.time()
  
  ## define parameters
  # fractional integration
  fractional_d <- 0 # in the open interval (-0.5, 0.5), positive number means long-range dependence
  
  # autoregressive operators
  operators_ar <- list(
    function(x,y){ 0.3*sin(x-y) },
    function(x,y){ 0.3*cos(x-y) },
    function(x,y){ 0.3*sin(2*x) },
    function(x,y){ 0.3*cos(y) }
  )
  
  # moving average kernels
  operators_ma <- list(
    function(x,y){ x+y },
    function(x,y){ x },
    function(x,y){ y }
  )
  
  # covariance of the inovation
  if (method_i <= 2){
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
    farima_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma,sigma_eigenvalues=sigma_eigenvalues,sigma_eigenfunctions=sigma_eigenfunctions)
  } else { # automatic
    
      sigma <- function(x,y){
        1*sin(2*pi*x)*sin(2*pi*y)+
          0.6*cos(2*pi*x)*cos(2*pi*y)+
          0.3*sin(4*pi*x)*sin(4*pi*y)+
          0.1*cos(4*pi*x)*cos(4*pi*y)+
          0.1*sin(6*pi*x)*sin(6*pi*y)+
          0.1*cos(6*pi*x)*cos(6*pi*y)+
          0.05*sin(8*pi*x)*sin(8*pi*y)+
          0.05*cos(8*pi*x)*cos(8*pi*y)+
          0.05*sin(10*pi*x)*sin(10*pi*y)+
          0.05*cos(10*pi*x)*cos(10*pi*y)
      }
    
    farima_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma,sigma=sigma)
  }
  
  
  # switch if to use fully spectral simulation or hybrid
  if (method_i %in% c(1,3)){ hybrid_ar <- FALSE }
  if (method_i %in% c(2,4)){ hybrid_ar <- TRUE }
  
  # simulate trajectory HERE
  fts_x <- FARFIMA_simulate(farima_pars, t_max, n_grid, seed_number=seed_number, hybrid_ar=hybrid_ar)
  
  # ending time
  end_time  = Sys.time()
  timing <- difftime(end_time,start_time, units="secs")  # saving the simulation time
  
  # evaluate and save empirical cov kernels
  for (lag_i in 1:lag_to_compare_n){
    lag <- lag_to_compare[lag_i]
    
    # calculate covariance kernel
    m <- cov( t(fts_x[,(lag+1):t_max]), t(fts_x[,1:(t_max-lag)]) )
    covlags_all[lag_i,,] <- m
  }
  
  return(list(covlags_all=covlags_all,timing=timing))
}


###################################################################################
## simulate ARMA(4,3)

simulation_run_ARMA_lowrank_space <- function( t_max, n_grid, lag_to_compare, seed_number=NaN ){
  
  # here I'm saving empirical col kernels for the designated lags, for each sim run
  lag_to_compare_n <- length(lag_to_compare)
  covlags_all <- array(NaN, dim=c(lag_to_compare_n,n_grid,n_grid) )
  
  # simulate trajectory
  start_time = Sys.time()
  
  ## define parameters
  # fractional integration
  fractional_d <- 0 # in the open interval (-0.5, 0.5), positive number means long-range dependence
  
  # autoregressive operators
  operators_ar <- list(
    function(x,y){ 0.3*sin(x-y) },
    function(x,y){ 0.3*cos(x-y) },
    function(x,y){ 0.3*sin(2*x) },
    function(x,y){ 0.3*cos(y) }
  )
  
  # moving average kernels
  operators_ma <- list(
    function(x,y){ x+y },
    function(x,y){ x },
    function(x,y){ y }
  )
  
  # covariance of the inovation
  sigma <- function(x,y){
    1*sin(2*pi*x)*sin(2*pi*y)+
      0.6*cos(2*pi*x)*cos(2*pi*y)+
      0.3*sin(4*pi*x)*sin(4*pi*y)+
      0.1*cos(4*pi*x)*cos(4*pi*y)+
      0.1*sin(6*pi*x)*sin(6*pi*y)+
      0.1*cos(6*pi*x)*cos(6*pi*y)+
      0.05*sin(8*pi*x)*sin(8*pi*y)+
      0.05*cos(8*pi*x)*cos(8*pi*y)+
      0.05*sin(10*pi*x)*sin(10*pi*y)+
      0.05*cos(10*pi*x)*cos(10*pi*y)
  }
  
  # put the parameters into one list
  farima_pars <- list(fractional_d=fractional_d,operators_ar=operators_ar,operators_ma=operators_ma,sigma=sigma)
  
  # simulate trajectory HERE
  fts_x <- simulate_ARMA_space(farima_pars,t_max,n_grid, seed_number = seed_number)
  
  # ending time
  end_time  = Sys.time()
  timing <- difftime(end_time,start_time, units="secs")  # saving the simulation time
  
  # evaluate and save empirical cov kernels
  for (lag_i in 1:lag_to_compare_n){
    lag <- lag_to_compare[lag_i]
    
    # calculate covariance kernel
    m <- cov( t(fts_x[,(lag+1):t_max]), t(fts_x[,1:(t_max-lag)]) )
    covlags_all[lag_i,,] <- m
  }
  
  return(list(covlags_all=covlags_all,timing=timing))
}


###################################################################################
## simulate FARIMA using Li et al.

simulation_run_FARIMA_Lietal <- function( t_max, n_grid, lag_to_compare, seed_number){
  
  # here I'm saving empirical col kernels for the designated lags, for each sim run
  lag_to_compare_n <- length(lag_to_compare)
  covlags_all <- array(NaN, dim=c(lag_to_compare_n,n_grid,n_grid) )
  
  # simulate trajectory
  start_time = Sys.time()
  fts_x <- sim_FARMA(n = t_max, d = 0.2, no_grid = n_grid, seed_number = seed_number, process="FAR")
  end_time = Sys.time()
  timing <- difftime(end_time,start_time, units="secs") # saving the simulation time
  
  # evaluate and save empirical cov kernels
  for (lag_i in 1:lag_to_compare_n){
    lag <- lag_to_compare[lag_i]
    
    # save
    m <- cov( t(fts_x[,(lag+1):t_max]), t(fts_x[,1:(t_max-lag)]) )
    #persp(m)
    covlags_all[lag_i,,] <- m
  }
  
  return(list(covlags_all=covlags_all,timing=timing))
}


###############################################################################################
calculate_errors_from_avg_cov <- function(covlags_all, n_grid, lag_to_compare, process){
  
  lag_to_compare_n <- length(lag_to_compare)
  
  # here I shall save my errors for each lag
  error_fro <- numeric(lag_to_compare_n)
  error_nuc <- numeric(lag_to_compare_n)
  error_sup <- numeric(lag_to_compare_n)
  
  # analyse my lags of interest
  for (lag_i in 1:lag_to_compare_n){
    lag <- lag_to_compare[lag_i]
    
    # true cov
    file_name <- paste("",process,"_covs/",process,"_lag_",lag,"_ngrid_",n_grid,".txt",sep="")
    covlag_true <- Re(as.matrix(read.table(file_name)))
    
    # average covariance for this lag
    covlag_average <- covlags_all[lag_i,,]
    
    # evaluate errors
    residuals <- covlag_true - covlag_average
    error_fro[lag_i] <- norm(residuals, type="f") / n_grid # frobenius norm
    error_nuc[lag_i] <- sum( svd(residuals,0,0)$d ) / n_grid
    error_sup[lag_i] <- max(abs(residuals)) # supremum norm
  }
  
  # lag-0 operator
  file_name <- paste("",process,"_covs/",process,"_lag_",0,"_ngrid_",n_grid,".txt",sep="")
  covlag0_true <- Re(as.matrix(read.table(file_name)))
  
  # relative errors w.r.t. the true lag-0 operator
  relative_error_fro <- error_fro / (norm(covlag0_true, type="f") / n_grid)
  relative_error_nuc <- error_nuc / ( sum(diag(covlag0_true)) / n_grid)
  relative_error_sup <- error_sup / max(abs(covlag0_true))
  
  return( list(error_fro=error_fro, relative_error_fro=relative_error_fro,
               error_nuc=error_nuc, relative_error_nuc=relative_error_nuc,
               error_sup=error_sup, relative_error_sup=relative_error_sup) )
}




###############################################################################################
simulate_ARMA_space <- function(arma_pars,t_max,n_grid, burn_in = 500, seed_number = NaN){
  
  ## random seed if assigned
  if (!is.nan(seed_number)){ set.seed(seed_number) }
  
  # grid for evaluation
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # save model order
  ar_order <- length(arma_pars$operators_ar)
  ma_order <- length(arma_pars$operators_ma)
  
  # firstly simply copy the data structure
  arma_pars_space <- arma_pars
  
  # express sigma
  arma_pars_space$sigma <- arma_pars$sigma(grid_matrix,t(grid_matrix))
  
  # express AR operators
  if (ar_order > 0){
    for (j in 1:ar_order){
      arma_pars_space$operators_ar[[j]] <- arma_pars$operators_ar[[j]](grid_matrix,t(grid_matrix))
    }
  }
  
  # express MA operators
  if (ma_order > 0){
    for (j in 1:ma_order){
      arma_pars_space$operators_ma[[j]] <- arma_pars$operators_ma[[j]](grid_matrix,t(grid_matrix))
    }
  }
  
  # the FTS data structure
  fts_x <- matrix(0, nrow = n_grid, ncol = t_max + burn_in)
  
  # define the innovation noise
  epsilon <- t(rmvnorm(t_max + burn_in, sigma = arma_pars_space$sigma))
  
  # AR cycle
  for (t in (max(ma_order,ar_order)+1):(t_max + burn_in)){
    
    # AR part
    if (ar_order > 0){
      for (ii in 1:ar_order){
        fts_x[,t] <- fts_x[,t] + t(arma_pars_space$operators_ar[[ii]] %*% fts_x[,t-ii]) / n_grid
      }
    }
    
    # MA part
    fts_x[,t] <- fts_x[,t] + epsilon[,t]
    if (ma_order > 0){
      for (ii in 1:ma_order){
        fts_x[,t] <- fts_x[,t] + t(arma_pars_space$operators_ma[[ii]] %*% epsilon[,t-ii]) / n_grid
      }
    }
  }
  
  return( fts_x[,(burn_in+1):(burn_in+t_max)] )
  
}

