

calculate_lagh <- function(n_grid, lag){
  ##################################################
  ## precalculation setting
  
  grid <- seq( 0, 1, length.out = n_grid ) # grid of [0,1] interval
  
  # parameters
  n_pc <- 1000
  fractional_d <- 0.2
  
  ##################################################
  ## edit R's function to allow for complex valued integration
  integrate_complex <- function(fun,...){
    fun_re <- function(x){ Re(fun(x)) }
    fun_im <- function(x){ Im(fun(x)) }
    cutoff <- 0.0001
    
    integral_re <- integral(fun_re, xmin=0, xmax=cutoff) +
      integral(fun_re, xmin=cutoff, xmax=2*pi-cutoff) +
      integral(fun_re, xmin=2*pi-cutoff, xmax=2*pi)
    
    return( integral_re )
  }
  
  
  ##################################################
  ## start integration
  
  start_time = Sys.time()
  # where shall I save the cov lag h kernel
  covlagh = matrix(0, nrow=n_grid, ncol=n_grid)
  
  # prepare quantities for calculation of basis functions
  exp_f <- exp( grid^2 / 2 )
  exp_f_norm2 <- sqrt(pi)/2 * erfi(1) # the norm of exp_f (calculated analytically)
  
  # integrate the spectral density
  pb <- txtProgressBar(style = 3)
  for (n in (1:n_pc)){
    setTxtProgressBar(pb, n/n_pc)
    # prepare eigenvalue and eigenfunction of Sigma
    sigma_l_sqrt <- 1/((n-0.5)*pi)
    sigma_e <- sqrt(2) * sin( (n-0.5)*pi*grid )
    
    f_inner_product <- function(x) { exp((x^2)/2) * sqrt(2) * sin( (n-0.5)*pi*x ) }
    inner_product_e_f <- integrate(f_inner_product, lower = 0, upper = 1, subdivisions=2000)$value
    #inner_product_e_f <- (exp(((1 - 2*n)^2 *pi^2)/8)*sqrt(pi)*
    #                        (-2*erfz(((-1 + 2*n)*pi)/(2*sqrt(2))) + erfz((-2i - pi + 2*n*pi)/(2*sqrt(2))) + erfz((2i - pi + 2*n*pi)/(2*sqrt(2)))))/2;
    #(E^(((1 - 2 n)^2 Pi^2)/8) Sqrt[Pi] (-2 Erf[((-1 + 2 n) Pi)/(2 Sqrt[2])] + Erf[(-2 I - Pi + 2 n Pi)/(2 Sqrt[2])] + Erf[(2 I - Pi + 2 n Pi)/(2 Sqrt[2])]))/2
    
    # define functions that I'm going to numerically integrate
    # c_omega <- function(omega){ ( exp(-1i*omega) *0.34 )/( 1- exp(-1i*omega)*0.34*exp_f_norm2 ) }
    c_omega <- function(omega){ ( exp(-1i*omega)*0.34 - 0.34^2 * exp_f_norm2  )/( 1 + (0.34*exp_f_norm2)^2 - 2*cos(omega)*0.34*exp_f_norm2 ) }
    f_integral_1 <- function(omega) { sigma_l_sqrt^2 * exp(1i*omega*lag) * (1/(2*pi)) * ( 2*sin(omega/2) )^(-2*fractional_d) }
    f_integral_2 <- function(omega) { sigma_l_sqrt^2 * inner_product_e_f * c_omega(omega) * exp(1i*omega*lag)*(1/(2*pi)) * ( 2*sin(omega/2) )^(-2*fractional_d)}
    f_integral_3 <- function(omega) { sigma_l_sqrt^2 * Conj(inner_product_e_f) * Conj(c_omega(omega)) * exp(1i*omega*lag)*(1/(2*pi)) * ( 2*sin(omega/2) )^(-2*fractional_d) }
    f_integral_4 <- function(omega) { sigma_l_sqrt^2 * abs(inner_product_e_f)^2 * abs(c_omega(omega))^2 * exp(1i*omega*lag)*(1/(2*pi)) * ( 2*sin(omega/2) )^(-2*fractional_d) }
    
    integral_1 <- integrate_complex(f_integral_1)
    integral_2 <- integrate_complex(f_integral_2)
    integral_3 <- integrate_complex(f_integral_3)
    integral_4 <- integrate_complex(f_integral_4)
    
    
    covlagh <- covlagh + sigma_e %o% (sigma_e) *integral_1
    covlagh <- covlagh + exp_f %o% (sigma_e) * integral_2
    covlagh <- covlagh + sigma_e %o% (exp_f) * integral_3
    covlagh <- covlagh + exp_f %o% (exp_f) * integral_4
    
  }
  close(pb)
  
  return(covlagh)
}


for (n_grid in c(101)){
  for (lag in c(0,1,2,3,5,10,20,30,40,60,80,100)){
    covlagh <- calculate_lagh(n_grid, lag)
    name <- paste("farima_covs/farima_lag_",lag,"_ngrid_",n_grid,".txt",sep="")
    print(name)
    write.table(Re(covlagh), file = name)
  }
}

for (n_grid in c(101,201,501,1001)){
  for (lag in c(0)){
    covlagh <- calculate_lagh(n_grid, lag)
    name <- paste("farima_covs/farima_lag_",lag,"_ngrid_",n_grid,".txt",sep="")
    print(name)
    write.table(Re(covlagh), file = name)
  }
}

