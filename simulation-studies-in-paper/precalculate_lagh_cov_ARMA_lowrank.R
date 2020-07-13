setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")


library(sde)
library(pracma)
library(mvtnorm)




f_calculate_lagh <- function(n_grid, lags_all){
  ##################################################
  ## precalculation setting
  
  # frequency intergration grid 
  n_grid_freq <- 2000
  grid_freq <- seq( 0, pi, length.out = n_grid_freq )
  
  
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
  
  
  
  lags_all_n <- length(lags_all)

  ######################################
  # evaluate operators on grid
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  
  # save model order
  ar_order <- length(operators_ar)
  ma_order <- length(operators_ma)
  
  # express sigma
  sigma_eval <- sigma(grid_matrix,t(grid_matrix))  / n_grid
  
  # express AR operators
  operators_ar_eval <- operators_ar
  if (ar_order > 0){
    for (j in 1:ar_order){
      operators_ar_eval[[j]] <- operators_ar[[j]](grid_matrix,t(grid_matrix)) / n_grid
    }
  }
  
  # express MA operators
  operators_ma_eval <- operators_ma
  if (ma_order > 0){
    for (j in 1:ma_order){
      operators_ma_eval[[j]] <- operators_ma[[j]](grid_matrix,t(grid_matrix))  / n_grid
    }
  }
  
  ##################################################
  ## start integration
  cov_lags <- array(0, dim=c(lags_all_n,n_grid,n_grid) )
  
  pb <- txtProgressBar(style = 3)
  for (k in 1:n_grid_freq){
  # cov_lags<-foreach(k = 1:n_grid_freq, .combine='+', .packages=c('pracma','sde','mvtnorm')) %dopar% {
    setTxtProgressBar(pb, k/n_grid_freq)
    omega <- grid_freq[k]
    
    # evaluate spectral density
    # express AR operators
    ar_part <- diag(n_grid)
    if (ar_order > 0){
      for (j in 1:ar_order){
        ar_part <- ar_part - operators_ar_eval[[j]] * exp(-1i*j*omega)
      }
    }
    ar_part_inv <- solve(ar_part)
    
    # express MA operators
    ma_part <- diag(n_grid)
    if (ma_order > 0){
      for (j in 1:ma_order){
        ma_part <- ma_part + operators_ma_eval[[j]] * exp(-1i*j*omega)
      }
    }
    
    # spec_density
    spec_density <- 1/(2*pi) * ar_part_inv %*% ma_part %*% sigma_eval %*% Conj(t(ma_part)) %*% Conj(t(ar_part_inv)) * n_grid
    
    # contributions to cov lags
    for (lag_i in 1:lags_all_n){
      lag <- lags_all[lag_i]
      cov_lags[lag_i,,] <- cov_lags[lag_i,,] + pi / n_grid_freq * (spec_density * exp(1i*lag*omega) + t(spec_density) * exp(-1i*lag*omega) )
    }
    
    
  }
  # close(pb)
  
  return(cov_lags)
  
}

for (n_grid in c(101,201,501,1001)){
  
  if (n_grid == 101){
    lags_all <- c(0,1,2,3,5,10,20,30,40,60,80,100)
  } else {
    lags_all <- c(0)
  }

  cov_lags <- f_calculate_lagh(n_grid, lags_all)

  # save lags
  for (lag_i in 1:length(lags_all)){
    lag <- lags_all[lag_i]
    name <- paste("arma_lowrank_covs/arma_lowrank_lag_",lag,"_ngrid_",n_grid,".txt",sep="")
    print(name)
    write.table(Re(cov_lags[lag_i,,]), file = name)
  }

}


