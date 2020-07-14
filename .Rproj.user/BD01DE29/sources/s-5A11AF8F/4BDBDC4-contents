
############################################
## packages loading and setwd

setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")
source("functions_simstudy.R")

library(sde)
library(pracma)
library(mvtnorm)
library(specsimfts)

# global settings
# global_setting_cores <- 2
# registerDoParallel(global_setting_cores)


# distribute all runs
t_max_all <- c()
n_grid_all <- c()
lag_to_compare_all <- list()
method_i_all <- c()
file_name_all <- c()
sim_runs <- 0



n_pc <- 101
t_max_all <- c(50,100,200,400,800,1600,3200,6400)
n_grid <- 101


k_bbridge <- function(x,y) { pmin(x,y)-x*y }
spec_density <- function( omega, x,y ){ 1/(1-0.9*cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }

harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*((x-omega/pi)%%1 ))  ) }


#################################################################################################################################
# using CKL
timing_CKL <- NULL
for (t_max_i in 1:length(t_max_all)){
  t_max <- t_max_all[t_max_i]
  
  start_time = Sys.time()
  fts_x <- HKL_simulate(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc)
  end_time  = Sys.time()
  timing_CKL[t_max_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}


#################################################################################################################################
# using numerical eigen-decomposition
timing_SVD <- NULL
for (t_max_i in 1:length(t_max_all)){
  t_max <- t_max_all[t_max_i]
  
  start_time = Sys.time()
  n_pc <- n_grid
  fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc)
  end_time  = Sys.time()
  timing_SVD[t_max_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}


print(timing_CKL)
print(timing_SVD)











