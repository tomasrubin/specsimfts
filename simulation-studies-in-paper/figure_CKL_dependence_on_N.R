setwd("C:/Users/tomas/OneDrive/Documents/GitHub/specsimfts/simulation-studies-in-paper/")

library(specsimfts)
library(RColorBrewer)
library(pracma)



harmonic_eigenvalues <- function( omega, n ){ 1/( (1-0.9 *cos(omega)) * (n*pi)^2 ) }
harmonic_eigenfunctions <- function(omega, n, x){ sqrt(2)*sin( n*(pi*x-omega)  ) }


## simulation setting
seed <- 5
t_max <- 100
n_grid <- 1001
n_pc_all <- c(1,2,3,5,10,25,50,1001)
fts_x_all <- matrix(NA, ncol=length(n_pc_all), nrow=n_grid)
fts_x_profile <- matrix(NA, ncol=length(n_pc_all), nrow=100)

# simulate trajectories
for (n_pc_i in 1:length(n_pc_all)){
  n_pc <- n_pc_all[n_pc_i]
  fts_x <- HKL_simulate(harmonic_eigenvalues, harmonic_eigenfunctions, t_max, n_grid, n_pc, seed_number = seed, include_freq_zero = T)
  fts_x_all[,n_pc_i] <- fts_x[,1]
  fts_x_profile[,n_pc_i] <- fts_x[ceil(n_grid/2),1:100]
}


cols <- brewer.pal(8, "Dark2")

plot( seq(0,1,length.out=n_grid), fts_x_all[,length(n_pc_all)], type="l", col=cols[length(n_pc_all)], lwd=2, main="custom defined harmonic Karhunen-Loeve expansion", ylab="",xlab="") # , ylim=c(-0.65,0.55)
for (n_pc_i in rev(1:(length(n_pc_all)-1))){
  lines(seq(0,1,length.out=n_grid),fts_x_all[,n_pc_i], col=cols[n_pc_i], lwd=2)
}


legend( x="top", legend=n_pc_all, col=cols, lty = 1, lwd=2, ncol=8 )
