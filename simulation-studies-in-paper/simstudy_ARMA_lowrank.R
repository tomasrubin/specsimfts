#################################################################################
## simulation study for the lowrank ARMA process considered in the paper
## running of this script takes hours if not a day
#################################################################################

setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")
source("functions_simstudy.R")


library(sde)
library(pracma)
library(foreach)
library(doParallel)
library(mvtnorm)
library(specsimfts)


# global settings
global_setting_sim_run_n <- 100
global_setting_cores <- 7
registerDoParallel(global_setting_cores)
methods <- c("lowrank_spec","lowrank_hybrid","svd_spec","svd_hybrid","space")
methods_run_now <- 1:5


# distribute all runs
t_max_all <- c()
n_grid_all <- c()
lag_to_compare_all <- list()
method_i_all <- c()
file_name_all <- c()
sim_runs <- 0



# investigate n_grid
for (n_grid in rev(c(101,201,501,1001))){# (c(101,201,501,1001,2001))){
  for (ii in methods_run_now){
    sim_runs <- sim_runs+1
    t_max_all[sim_runs] <- 800
    n_grid_all[sim_runs] <- n_grid
    lag_to_compare_all[[sim_runs]] <- c(0)
    method_i_all[sim_runs] <- ii
    file_name_all[sim_runs] <- "results/ARMA_lowrank_n_grid.csv"
  }
}

# t_max investigate
for (t in rev(c(50,100,200,400,800,1600,3200,6400))){
  for (ii in methods_run_now){
    sim_runs <- sim_runs+1
    t_max_all[sim_runs] <- t
    n_grid_all[sim_runs] <- 101
    lag_to_compare_all[[sim_runs]] <- c(0)
    method_i_all[sim_runs] <- ii
    file_name_all[sim_runs] <- "results/ARMA_lowrank_t_max.csv"
  }
}

# lags
for (ii in methods_run_now){
  sim_runs <- sim_runs+1
  t_max_all[sim_runs] <- 800
  n_grid_all[sim_runs] <- 101
  lag_to_compare_all[[sim_runs]] <- c(0,1,2,3,5,10,20,30,40,60,80,100)
  method_i_all[sim_runs] <- ii
  file_name_all[sim_runs] <- "results/ARMA_lowrank_lags.csv"
}



############################################################################################################################################
## precalculate autocovariance operators
source("precalculate_lagh_cov_ARMA_lowrank.R")

############################################################################################################################################
## run it

foreach(run_i = 1:sim_runs, .combine=c, .packages=c('pracma','sde','mvtnorm')) %dopar% {
# for(run_i in 1:sim_runs){

  
  t_max <- t_max_all[run_i]
  n_grid <- n_grid_all[run_i]
  method_i <- method_i_all[run_i]
  lag_to_compare <- lag_to_compare_all[[run_i]]
  file_name <- file_name_all[run_i]
  
  # structure to save results
  covlags_avg <- array(0, dim=c(length(lag_to_compare),n_grid,n_grid) )
  timing_avg <- 0
  
  # cycle for simulation runs
  for (sim_run_i in 1:global_setting_sim_run_n){
    # print(sim_run_i)
    
    seed_number <- run_i * 1000 + sim_run_i
    
    if (method_i <= 4){
      run <- simulation_run_ARMA_lowrank( t_max, n_grid, lag_to_compare, method_i=method_i, seed_number=seed_number )
    } else {
      run <- simulation_run_ARMA_lowrank_space( t_max, n_grid, lag_to_compare, seed_number=seed_number )
    }
    
    covlags_avg <- covlags_avg + run$covlags_all / global_setting_sim_run_n
    timing_avg <- timing_avg + run$timing / global_setting_sim_run_n
  }
  
  
  # get simulation error by average empirical covariance
  errors <- calculate_errors_from_avg_cov(covlags_avg, n_grid, lag_to_compare, process="arma_lowrank")
  
  ## save resultsprepare date frame to save
  dto_save <- data.frame( method=methods[method_i], n_simul=global_setting_sim_run_n, t_max, n_grid, time=as.numeric(timing_avg), lag=lag_to_compare,
                            error_fro = errors$error_fro, rel_error_fro = errors$relative_error_fro,
                            error_nuc = errors$error_nuc, rel_error_nuc =  errors$relative_error_nuc,
                            error_sup = errors$error_sup, rel_error_sup =  errors$relative_error_sup)
  
  
  if (file.exists(file_name)){
    write.table(dto_save, file=file_name, sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)  
  } else {
    write.table(dto_save, file=file_name, sep = ",", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  
}













