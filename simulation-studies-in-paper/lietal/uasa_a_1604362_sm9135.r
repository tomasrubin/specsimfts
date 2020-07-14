##############################
# install and load R packages
##############################

install.packages(c("sde", "doMC", "sandwich", "ggplot2", "gridExtra"))

require(sde)
require(doMC)
require(sandwich)
require(ggplot2)
require(gridExtra)

#############################
# set up a working directory
#############################

setwd("simulation")
source("sim_FARMA.R")

# estimate H value, d = H - 0.5

H_estimate <- function(fun_dat, band_const = 1, choice_score = c("first", "average"), 
                       est_method = c("classical", "modified"))
{
    est_method = match.arg(est_method)
    T = ncol(fun_dat)
    m = min(T - 1, nrow(fun_dat)) * band_const

    X_bar = rowMeans(fun_dat, na.rm = TRUE)
    center_dat = sweep(fun_dat, 1, X_bar)
   
    # calculating long-run covariance for a given lag
    
    gamma_l <- function(lag, T)
    {
        gamma_lag_sum = 0
        if(lag >= 0)
        {
            for(ij in 1:(T-lag))
            {
                gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij+lag)])))
            }
        }
        else
        {
            for(ij in 1:(T - abs(lag)))
            {
                gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij + abs(lag)]) %*%
                                t(as.matrix(center_dat[,ij])))
            }
        }
        return(gamma_lag_sum/T)
    }
   
    # calculating the sum of long-run covariance for all lags (m equals the total number of grid points)
    
    C = 0
    index = seq(-m, m, by = 1)
    for(k in 1:length(index))
    {
        C = C + (m - abs(index[k])) * gamma_l(lag = index[k], T = T)
    }
    
    # eigen decomposition
    
    eigen_decomp = eigen(C)   
    
    # select the number of component based on 95% total amount of variation
    
    prop_eigen_value = cumsum(eigen_decomp$values/sum(eigen_decomp$values))
    ncomp = head(which(prop_eigen_value >= 0.95), 1)
   
    # first functional principal component
    
    if(choice_score == "average")
    {
        if(ncomp == 1)
        {
            eigen_function = matrix(eigen_decomp$vectors[,1], nrow = 1)
            score = as.numeric(eigen_function %*% fun_dat)
        }
        if(ncomp > 1)
        {
            eigen_function = eigen_decomp$vectors[,1:ncomp]
        
            # averaging principal component scores
            score = colMeans(t(eigen_function) %*% fun_dat)
        }
    }
    if(choice_score == "first")
    {
        eigen_function = matrix(eigen_decomp$vectors[,1], nrow = 1)
        score = as.numeric(eigen_function %*% fun_dat)
    }
    
    # calculating R_n and S_n values
    
    R_value = var_value = vector("numeric",length(score))
    for(ik in 1:length(score))
    {
        R_value[ik] = sum(score[1:ik] - mean(score))
        var_value[ik] = (score[ik] - mean(score))^2
    }
    R_n = max(R_value) - min(R_value)
    
    S_n_classic = sqrt(sum(var_value)/length(score))
    if(est_method == "classical")
    {
        R_over_S = R_n/S_n_classic
        tilde_H1 = log(R_over_S)/log(length(score))
        d_est = tilde_H1 - 0.5
    }
    if(est_method == "modified")
    {
        rho_hat = acf(score, type = "correlation", plot = FALSE)$acf[2,,1]
        q = floor((length(score) * 3 / 2)^(1/3) * ((2 * rho_hat/(1 - rho_hat^2))^(2/3)))
    
        rho_cov = acf(score, type = "covariance", plot = FALSE)$acf[2:(q+1),,1]
        weight = weight_multi = vector("numeric", q)
        for(ik in 1:q)
        {
            weight[ik] = 1 - ik/(q + 1)
            weight_multi[ik] = weight[ik] * rho_cov[ik] 
        }
        S_n = sqrt(sum(var_value)/length(score) + 2 * sum(weight_multi))
        tilde_H1 = R_n/(S_n * sqrt(length(score))) # for hypothesis test only
        d_est = tilde_H1 - 0.5
    }    
    
    # alpha_estimate
    
    alpha = 1.5 - tilde_H1
    return(list(C = C, score = score, H_est = tilde_H1, d_est = d_est, alpha_est = alpha, 
                ncomp = ncomp, prop_eigen_value = prop_eigen_value, 
                eigen_function = as.numeric(eigen_function)))
}

# a set of equally-spaced grid points bounded between 0 and 1

grid = seq(0, 1, length.out = 101)

# R = 1000 replications with different psuedo random seeds

d_fun <- function(ik, d_val, sample_size, score_choice, estimate_method, const_band)
{
    res = H_estimate(fun_dat = sim_FARMA(n = sample_size, d = d_val, no_grid = 101, seed_number = ik), 
                     band_const = const_band, choice_score = score_choice, est_method = estimate_method)
    return(c(res$d_est, res$ncomp))
}

d_fun_FARMA <- function(ik, d_val, sample_size, score_choice, estimate_method, const_band)
{
    res = H_estimate(fun_dat = sim_FARMA(n = sample_size, d = d_val, no_grid = 101, seed_number = ik, process = "FARMA"), 
                     band_const = const_band, choice_score = score_choice, est_method = estimate_method)
    return(c(res$d_est, res$ncomp))
}

##############
# FAR example
##############

# d = 0.2 (the one in the simulation study)

d_0.2_H_est_n500_example = H_estimate(fun_dat = sim_FARMA(n = 500, d = 0.2, no_grid = 101, seed_number = 1), choice_score = "first")
d_0.2_H_est_n1000_example = H_estimate(fun_dat = sim_FARMA(n = 1000, d = 0.2, no_grid = 101, seed_number = 1), choice_score = "first")
d_0.2_H_est_n2000_example = H_estimate(fun_dat = sim_FARMA(n = 2000, d = 0.2, no_grid = 101, seed_number = 1), choice_score = "first")

d_0.2_n500_example_plot1 = qplot(grid, as.numeric(d_0.2_H_est_n500_example$eigen_function),
                                 xlab = "Grid points", ylab = "First eigenfunction",
                                 main = "n = 500", ylim = c(-0.17,-0.01), geom = "line") 

d_0.2_n1000_example_plot2 = qplot(grid, as.numeric(d_0.2_H_est_n1000_example$eigen_function),
                                  xlab = "Grid points", ylab = "First eigenfunction",
                                  main = "n = 1000", ylim = c(-0.17,-0.01), geom = "line") 

d_0.2_n2000_example_plot3 = qplot(grid, as.numeric(d_0.2_H_est_n2000_example$eigen_function),
                                  xlab = "Grid points", ylab = "First eigenfunction",
                                  main = "n = 2000", ylim = c(-0.17,-0.01), geom = "line") 
grid.arrange(d_0.2_n500_example_plot1, d_0.2_n1000_example_plot2, d_0.2_n2000_example_plot3, ncol=3)

################
# FARMA example
################

d_0.2_H_est_n500_FARMA_example = H_estimate(fun_dat = sim_FARMA(n = 500, d = 0.2, no_grid = 101, seed_number = 1, process = "FARMA"), choice_score = "first")
d_0.2_H_est_n1000_FARMA_example = H_estimate(fun_dat = sim_FARMA(n = 1000, d = 0.2, no_grid = 101, seed_number = 1, process = "FARMA"), choice_score = "first")
d_0.2_H_est_n2000_FARMA_example = H_estimate(fun_dat = sim_FARMA(n = 2000, d = 0.2, no_grid = 101, seed_number = 1, process = "FARMA"), choice_score = "first")

d_0.2_FARMA_example_plot1 = qplot(grid, as.numeric(d_0.2_H_est_n500_FARMA_example$eigen_function),
                                  xlab = "Grid points", ylab = "First eigenfunction",
                                  main = "n = 500", ylim = c(-0.16,-0.02), geom = "line")

d_0.2_FARMA_example_plot2 = qplot(grid, as.numeric(d_0.2_H_est_n1000_FARMA_example$eigen_function),
                                  xlab = "Grid points", ylab ="First eigenfunction",
                                  main = "n = 1000", ylim = c(-0.16,-0.02), geom = "line")

d_0.2_FARMA_example_plot3 = qplot(grid, as.numeric(d_0.2_H_est_n2000_FARMA_example$eigen_function),
                                  xlab = "Grid points", ylab = "First eigenfunction",
                                  main = "n = 2000", ylim = c(-0.16,-0.02), geom = "line")

grid.arrange(d_0.2_FARMA_example_plot1, d_0.2_FARMA_example_plot2, d_0.2_FARMA_example_plot3, ncol=3)

############################################################
## Monte-Carlo simulation for 1000 random samples (d = 0.2)
############################################################

## FAR (const_band = 1)

# set up the number of cores

registerDoMC(20)
first_d_0.2_H_est_n500 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 500, 
											estimate_method = "classical", score_choice = "first", const_band = 1)

registerDoMC(20)
first_d_0.2_H_est_n1000 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 1000, 
											estimate_method = "classical", score_choice = "first", const_band = 1)

registerDoMC(20)
first_d_0.2_H_est_n2000 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 2000, 
											estimate_method = "classical", score_choice = "first", const_band = 1)

## FARMA (const_band = 1)

registerDoMC(20)
first_d_0.2_H_est_n500_FARMA = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 500, score_choice = "first")


registerDoMC(20)
first_d_0.2_H_est_n1000_FARMA = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 1000, score_choice = "first")


registerDoMC(20)
first_d_0.2_H_est_n2000_FARMA = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 2000, score_choice = "first")


##########################################################
# sensitivity analysis for different bandwidth parameters
##########################################################

## FAR (const_band = 2)

registerDoMC(20)
first_d_0.2_H_est_n500_bandconst_2 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 500, score_choice = "first", 
                                                                         estimate_method = "classical", const_band = 2)

registerDoMC(20)
first_d_0.2_H_est_n1000_bandconst_2 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 1000, score_choice = "first", 
                                                                         estimate_method = "classical", const_band = 2)

registerDoMC(20)
first_d_0.2_H_est_n2000_bandconst_2 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 2000, score_choice = "first", 
                                                                         estimate_method = "classical", const_band = 2)

# FARMA (const_band = 2)

registerDoMC(20)
first_d_0.2_H_est_n500_FARMA_bandconst_2 = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 500, score_choice = "first",
                                                                                     estimate_method = "classical", const_band = 2)

registerDoMC(20)
first_d_0.2_H_est_n1000_FARMA_bandconst_2 = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 1000, score_choice = "first",
                                                                                      estimate_method = "classical", const_band = 2)

registerDoMC(20)
first_d_0.2_H_est_n2000_FARMA_bandconst_2 = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 2000, score_choice = "first",
                                                                                      estimate_method = "classical", const_band = 2)

##########################################################
# sensitivity analysis for different bandwidth parameters
##########################################################

# FAR (const_band = 0.5)

registerDoMC(20)
first_d_0.2_H_est_n500_bandconst_0.5 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 500, score_choice = "first", 
                                                                         estimate_method = "classical", const_band = 0.5)

registerDoMC(20)
first_d_0.2_H_est_n1000_bandconst_0.5 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 1000, score_choice = "first", 
                                                                          estimate_method = "classical", const_band = 0.5)

registerDoMC(20)
first_d_0.2_H_est_n2000_bandconst_0.5 = foreach(iwk = 1:1000) %dopar% d_fun(ik = iwk, d_val = 0.2, sample_size = 2000, score_choice = "first", 
                                                                          estimate_method = "classical", const_band = 0.5)

# FARMA (const_band = 0.5)

registerDoMC(20)
first_d_0.2_H_est_n500_FARMA_bandconst_0.5 = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 500, score_choice = "first",
                                                                                     estimate_method = "classical", const_band = 0.5)

registerDoMC(20)
first_d_0.2_H_est_n1000_FARMA_bandconst_0.5 = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 1000, score_choice = "first",
                                                                                      estimate_method = "classical", const_band = 0.5)

registerDoMC(20)
first_d_0.2_H_est_n2000_FARMA_bandconst_0.5 = foreach(iwk = 1:1000) %dopar% d_fun_FARMA(ik = iwk, d_val = 0.2, sample_size = 2000, score_choice = "first",
                                                                                      estimate_method = "classical", const_band = 0.5)

##########################
# summarizing the results
##########################

# FAR

# bandconst = 1

first_d_0.2_H_est_mat = cbind(do.call(rbind, first_d_0.2_H_est_n500)[,1],
                              do.call(rbind, first_d_0.2_H_est_n1000)[,1],
                              do.call(rbind, first_d_0.2_H_est_n2000)[,1])
colnames(first_d_0.2_H_est_mat) = c("n = 500","n = 1000","n = 2000")

ncomp_d_0.2_H_est_mat = cbind(do.call(rbind, first_d_0.2_H_est_n500)[,2],
                              do.call(rbind, first_d_0.2_H_est_n1000)[,2],
                              do.call(rbind, first_d_0.2_H_est_n2000)[,2])
colnames(ncomp_d_0.2_H_est_mat) = c("n = 500", "n = 1000", "n = 2000")

# bandconst = 2

first_d_0.2_H_est_mat_bandconst_2 = cbind(do.call(rbind, first_d_0.2_H_est_n500_bandconst_2)[,1],
                                          do.call(rbind, first_d_0.2_H_est_n1000_bandconst_2)[,1],
                                          do.call(rbind, first_d_0.2_H_est_n2000_bandconst_2)[,1])
colnames(first_d_0.2_H_est_mat_bandconst_2) = c("n = 500","n = 1000","n = 2000")

# bandconst = 0.5

first_d_0.2_H_est_mat_bandconst_0.5 = cbind(do.call(rbind, first_d_0.2_H_est_n500_bandconst_0.5)[,1],
                                            do.call(rbind, first_d_0.2_H_est_n1000_bandconst_0.5)[,1],
                                            do.call(rbind, first_d_0.2_H_est_n2000_bandconst_0.5)[,1])
colnames(first_d_0.2_H_est_mat_bandconst_0.5) = c("n = 500","n = 1000","n = 2000")


# FARMA

# bandconst = 1

first_d_0.2_H_est_mat_FARMA = cbind(do.call(rbind,first_d_0.2_H_est_n500_FARMA)[,1],
                                    do.call(rbind,first_d_0.2_H_est_n1000_FARMA)[,1],
                                    do.call(rbind,first_d_0.2_H_est_n2000_FARMA)[,1])
colnames(first_d_0.2_H_est_mat_FARMA) = c("n = 500", "n = 1000", "n = 2000")

ncomp_d_0.2_H_est_mat_FARMA = cbind(do.call(rbind,first_d_0.2_H_est_n500_FARMA)[,2],
                                    do.call(rbind,first_d_0.2_H_est_n1000_FARMA)[,2],
                                    do.call(rbind,first_d_0.2_H_est_n2000_FARMA)[,2])
colnames(ncomp_d_0.2_H_est_mat_FARMA) = c("n = 500", "n = 1000", "n = 2000")

# bandconst = 2

first_d_0.2_H_est_mat_FARMA_bandconst_2 = cbind(do.call(rbind, first_d_0.2_H_est_n500_FARMA_bandconst_2)[,1],
                                                do.call(rbind, first_d_0.2_H_est_n1000_FARMA_bandconst_2)[,1],
                                                do.call(rbind, first_d_0.2_H_est_n2000_FARMA_bandconst_2)[,1])
colnames(first_d_0.2_H_est_mat_FARMA_bandconst_2) = c("n = 500","n = 1000","n = 2000")

# bandconst = 0.5

first_d_0.2_H_est_mat_FARMA_bandconst_0.5 = cbind(do.call(rbind, first_d_0.2_H_est_n500_FARMA_bandconst_0.5)[,1],
                                                  do.call(rbind, first_d_0.2_H_est_n1000_FARMA_bandconst_0.5)[,1],
                                                  do.call(rbind, first_d_0.2_H_est_n2000_FARMA_bandconst_0.5)[,1])
colnames(first_d_0.2_H_est_mat_FARMA_bandconst_0.5) = c("n = 500","n = 1000","n = 2000")


round(rbind(colMeans(first_d_0.2_H_est_mat), apply(first_d_0.2_H_est_mat, 2, median)), 4)
round(rbind(colMeans(first_d_0.2_H_est_mat_bandconst_2), apply(first_d_0.2_H_est_mat_bandconst_2, 2, median)), 4)
round(rbind(colMeans(first_d_0.2_H_est_mat_bandconst_0.5), apply(first_d_0.2_H_est_mat_bandconst_0.5, 2, median)), 4)

round(rbind(colMeans(first_d_0.2_H_est_mat_FARMA), apply(first_d_0.2_H_est_mat_FARMA, 2, median)), 4)
round(rbind(colMeans(first_d_0.2_H_est_mat_FARMA_bandconst_2), apply(first_d_0.2_H_est_mat_FARMA_bandconst_2, 2, median)), 4)
round(rbind(colMeans(first_d_0.2_H_est_mat_FARMA_bandconst_0.5), apply(first_d_0.2_H_est_mat_FARMA_bandconst_0.5, 2, median)), 4)


# FAR

raw_data = as.numeric(as.matrix(first_d_0.2_H_est_mat))
sample_size = rep(c(rep(500, 1000), rep(1000, 1000), rep(2000, 1000)),3)
d_val = c(rep(0.1, 3000), rep(0.25, 3000), rep(0.4, 3000))

first_d_H_est_mat_d2 = data.frame(cbind(raw_data, sample_size, d_val))
first_d_H_est_mat_d2$sample = factor(first_d_H_est_mat_d2$sample_size, labels = c("500", "1000", "2000"))
first_d_H_est_mat_d2$d_val = factor(first_d_H_est_mat_d2$d_val, labels = c("d = 0.1", "d = 0.25", "d = 0.4"))

qplot(sample, raw_data, data = first_d_H_est_mat_d2, geom = "boxplot") + xlab("Sample size") + facet_grid(.~d_val,scales="free",space="free")+theme_bw() + ylab(expression(hat(d)))+ylim(c(-0.03,0.4))
boxplot(first_d_0.2_H_est_mat, ylim=c(0,0.3), ylab=expression(hat(d)))

# FARMA

raw_FARMA_data = as.numeric(as.matrix(first_d_H_est_mat_FARMA))

first_d_H_est_mat_d3 = data.frame(cbind(raw_FARMA_data, sample_size, d_val))
first_d_H_est_mat_d3$sample = factor(first_d_H_est_mat_d2$sample_size, labels = c("500", "1000", "2000"))
first_d_H_est_mat_d3$d_val = factor(first_d_H_est_mat_d2$d_val, labels = c("d = 0.1", "d = 0.25", "d = 0.4"))

qplot(sample, raw_FARMA_data, data = first_d_H_est_mat_d3, geom = "boxplot") + xlab("Sample size") + facet_grid(.~d_val,scales="free",space="free")+theme_bw() + ylab(expression(hat(d)))+ylim(c(0,0.4))
boxplot(first_d_0.2_H_est_mat_FARMA, ylim=c(0,0.3), ylab=expression(hat(d)))

# Compare FAR and FARMA models with d = 0.2 for different bandwidth constant values

par(mfrow=c(2,3))
boxplot(first_d_0.2_H_est_mat_bandconst_0.5, xlab="FAR model", ylab="d = 0.2")
title("0.5 x min(n-1, p)")
boxplot(first_d_0.2_H_est_mat, xlab="FAR model", ylab="d = 0.2")
title("min(n-1, p)")
boxplot(first_d_0.2_H_est_mat_bandconst_2, xlab="FAR model", ylab="d = 0.2")
title("2 x min(n-1, p)")

boxplot(first_d_0.2_H_est_mat_FARMA_bandconst_0.5, xlab="FARMA model", ylab="d = 0.2")
title("0.5 x min(n-1, p)")
boxplot(first_d_0.2_H_est_mat_FARMA, xlab="FARMA model", ylab="d = 0.2")
title("min(n-1, p)")
boxplot(first_d_0.2_H_est_mat_FARMA_bandconst_2, xlab="FARMA model", ylab="d = 0.2")
title("2 x min(n-1, p)")
