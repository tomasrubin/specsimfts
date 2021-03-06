#################################################################################
## simulation study for the lowrank ARMA process considered in the paper
## running of this script takes days
#################################################################################


############################################
## packages loading and setwd

setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")
source("functions_simstudy.R")

library(sde)
library(pracma)
library(foreach)
library(doParallel)
library(mvtnorm)
library(specsimfts)


# global settings
global_setting_sim_run_n <- 1000
global_setting_cores <- 4
registerDoParallel(global_setting_cores)




n_pc_all <- rev(c(1,2,3,5,10,20,30,50,100,200,300,500,1000))
t_max <- 1000
n_grid <- 1001
lag_to_compare <- c(0,1,2,3,5,10,20,30,40,60,80,100)

file_name <- "results/custom_CKL_1001.csv"

k_bbridge <- function(x,y) { pmin(x,y)-x*y }
spec_density <- function( omega, x,y ){ 1/(1-0.9*cos(omega)) * k_bbridge( (x-omega/pi)%%1, (y-omega/pi)%%1  ) }

harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*((x-omega/pi)%%1 ))  ) }



############################################################################################################################################
## precalculate autocovariance operators
source("precalculate_lagh_cov_custom_CKL.R")

############################################################################################################################################
## run it

foreach(n_pc_i = 1:length(n_pc_all), .combine=c, .packages=c('pracma','sde','mvtnorm','specsimfts')) %dopar% {

  n_pc <- n_pc_all[n_pc_i]

  # structure to save results
  covlags_avg <- array(0, dim=c(length(lag_to_compare),n_grid,n_grid) )
  timing_avg <- 0
  
  # cycle for simulation runs
  for (sim_run_i in 1:global_setting_sim_run_n){
    # print(sim_run_i)
    
    seed_number <- n_pc_i * 1000 + sim_run_i
    set.seed(seed_number)
    
    start_time = Sys.time()
    fts_x <- HKL_simulate(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc)
    end_time  = Sys.time()
    timing <- difftime(end_time,start_time, units="secs")  # saving the simulation time
    
    for (lag_i in 1:length(lag_to_compare)){
      lag <- lag_to_compare[lag_i]
      me <- cov( t(fts_x[,(1+lag):t_max]), t(fts_x[,1:(t_max-lag)]))
      covlags_avg[lag_i,,] <- covlags_avg[lag_i,,] + me / global_setting_sim_run_n
    }
    
    timing_avg <- timing_avg + timing / global_setting_sim_run_n  # saving the simulation time
  }
  
  
  # get simulation error by average empirical covariance
  errors <- calculate_errors_from_avg_cov(covlags_avg, n_grid, lag_to_compare, process="custom_CKL")
  
  ## save resultsprepare date frame to save
  df_to_save <- data.frame( method="CKL", n_simul=global_setting_sim_run_n, n_pc, t_max, n_grid, time=as.numeric(timing_avg), lag=lag_to_compare,
                            error_fro = errors$error_fro, rel_error_fro = errors$relative_error_fro,
                            error_nuc = errors$error_nuc, rel_error_nuc =  errors$relative_error_nuc,
                            error_sup = errors$error_sup, rel_error_sup =  errors$relative_error_sup)
  
  
  if (file.exists(file_name)){
    write.table(df_to_save, file=file_name, sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)  
  } else {
    write.table(df_to_save, file=file_name, sep = ",", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  
}













