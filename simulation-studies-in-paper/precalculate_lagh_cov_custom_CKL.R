setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")


library(sde)
library(pracma)
library(foreach)
library(doParallel)
library(mvtnorm)

# multi core setting
global_setting_cores <- 2
registerDoParallel(global_setting_cores)



k_bbridge <- function(x,y) { pmin(x,y)-x*y }
spec_density <- function( omega, x,y ){ 1/(1-0.9*cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }

harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*((x-omega/pi)%%1 ))  ) }




for(n_grid in c(1001)){
  
  
  lags_all <- c(0,1,2,3,5,10,20,30,40,60,80,100)
  

  # save lags
  foreach (lag_i = 1:length(lags_all), .combine=c, .packages=c('pracma','sde','mvtnorm')) %dopar% {
    
    lag <- lags_all[lag_i]
    cov_lag <- spec_density_covlagh_operator(spec_density, lag, n_grid)
    name <- paste("custom_CKL_covs/custom_CKL_lag_",lag,"_ngrid_",n_grid,".txt",sep="")
    print(name)
    write.table(Re(cov_lag), file = name)
  }

}



