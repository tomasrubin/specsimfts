
############################################
## packages loading and setwd

setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")
source("functions_simstudy.R")

library(sde)
library(pracma)
library(mvtnorm)
library(specsimfts)

# distribute all runs
t_max_all <- c()
n_grid_all <- c()
lag_to_compare_all <- list()
method_i_all <- c()
file_name_all <- c()
sim_runs <- 0



n_pc <- 100
t_max <- 1000
n_grid_all <- c(101,201,501,701,1001)


k_bbridge <- function(x,y) { pmin(x,y)-x*y }
spec_density <- function( omega, x,y ){ 1/(1-0.9*cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }

harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*((x-omega/pi)%%1 ))  ) }



#################################################################################################################################
# using CKL
timing_CKL <- NULL
for (n_grid_i in 1:length(n_grid_all)){
  n_grid <- n_grid_all[n_grid_i]
  
  start_time = Sys.time()
  fts_x <- HKL_simulate(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc)
  end_time  = Sys.time()
  timing_CKL[n_grid_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}

#################################################################################################################################
# using numerical eigen-decomposition
timing_SVD_full <- NULL
for (n_grid_i in 1:length(n_grid_all)){
  n_grid <- n_grid_all[n_grid_i]
  
  start_time = Sys.time()
  fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc=n_grid)
  end_time  = Sys.time()
  timing_SVD_full[n_grid_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}

#################################################################################################################################
# using numerical eigen-decomposition (PC = 5)
timing_SVD_5 <- NULL
for (n_grid_i in 1:length(n_grid_all)){
  n_grid <- n_grid_all[n_grid_i]
  
  start_time = Sys.time()
  fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc=5)
  end_time  = Sys.time()
  timing_SVD_5[n_grid_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}

#################################################################################################################################
# using numerical eigen-decomposition (PC = 10)
timing_SVD_10 <- NULL
for (n_grid_i in 1:length(n_grid_all)){
  n_grid <- n_grid_all[n_grid_i]
  
  start_time = Sys.time()
  fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc=10)
  end_time  = Sys.time()
  timing_SVD_10[n_grid_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}

#################################################################################################################################
# using numerical eigen-decomposition (PC = 50)
timing_SVD_50 <- NULL
for (n_grid_i in 1:length(n_grid_all)){
  n_grid <- n_grid_all[n_grid_i]
  
  start_time = Sys.time()
  fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc=50)
  end_time  = Sys.time()
  timing_SVD_50[n_grid_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}

#################################################################################################################################
# using numerical eigen-decomposition (PC = 100)
timing_SVD_100 <- NULL
for (n_grid_i in 1:length(n_grid_all)){
  n_grid <- n_grid_all[n_grid_i]
  
  start_time = Sys.time()
  fts_x <- spec_density_simulate(spec_density, t_max, n_grid, n_pc=100)
  end_time  = Sys.time()
  timing_SVD_100[n_grid_i] <- difftime(end_time,start_time, units="secs")  # saving the simulation time
}




print(timing_CKL)
print(timing_SVD_5)
print(timing_SVD_10)
print(timing_SVD_50)
print(timing_SVD_100)
print(timing_SVD_full)











