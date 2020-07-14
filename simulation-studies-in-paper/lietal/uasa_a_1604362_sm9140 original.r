#################################
# standard Brownian Motion [0,1]
#################################

BrownMat <- function(N, refinement)
{
    mat <- matrix(nrow = refinement, ncol = N)
    c <- 1
    while(c <= N)
    {
        vec <- BM(N = refinement-1)
        mat[,c] <- vec
        c <- c+1
    }
    return(mat)
}

##############################################
# functional kernel function in the FAR model
##############################################

funKernel = function(ref)
{
    Mat = matrix(nrow=ref, ncol=ref)
    for(i in 1:ref)
    {
        for(j in 1:ref)
        {
            Mat[i,j] = 0.34 * exp(0.5 * ((i/ref)^2 + (j/ref)^2))
        }
    }
    return(Mat)
}

funIntegral = function(ref, Mat, X)
{	
    Mat = Mat %*% X
    return(Mat/ref)
}

funARMat = function(refinement, eta_star_val)
{
    Mat = funKernel(refinement)
    res <- matrix(nrow = refinement, ncol = ncol(eta_star_val))
    res[,1] = eta_star_val[,1]
    for(ik in 2:(ncol(eta_star_val)))
    {
        res[,ik] = funIntegral(refinement, Mat, eta_star_val[,(ik-1)]) + eta_star_val[,ik]
    }
    return(res)
}

MAphi = function(ref, k)
{
    Mat = matrix(nrow = ref, ncol = ref)
    for(i in 1:ref)
    {
        for(j in 1:ref)
        {
            Mat[i,j] = min(i,j)/ref
        }
    }
    return(k * Mat)
}

funMAMat = function(refinement, eta_star_val)
{
    Mat = MAphi(refinement, 1.5)
    res <- matrix(nrow = refinement, ncol = ncol(eta_star_val)-1)
    for(ik in 2:(ncol(eta_star_val)))
    {
        res[,ik-1] = funIntegral(refinement, Mat, eta_star_val[,(ik-1)]) + eta_star_val[,ik]
    }
    return(res)
}

###############################################
# Simulated functional version of ARMA process
###############################################

sim_FARMA <- function(n = 500, d = 0.1, no_grid = 101, seed_number, process = c("FAR", "FARMA"))
{
    process = match.arg(process)
    set.seed(123 + seed_number)
  
    # simulate stochastic process realizations 
  
    step_1 = BrownMat(N = 2 * n + 100, refinement = no_grid)

    # calculating beta values
    
    beta_val = vector("numeric", n)
    for(i in 1:n)
    {
        beta_val[i] = exp(lgamma(i + d) - lgamma(i + 1) - lgamma(d))
    }
    
    # Step 2
    
    eta_star = matrix(NA, no_grid, (n+100))
    for(t_val in (n+1):(2*n+100))
    {
        inner_sum = matrix(NA, no_grid, n)
        for(ik in 1:n)
        {
            inner_sum[,ik] = (beta_val[ik] * step_1[,t_val - ik])
        }
        eta_star[,t_val - n] = step_1[,t_val] + rowSums(inner_sum)
    }

    # Step 3
    
    if(process == "FAR")
    {
	    sim_X = funARMat(refinement = no_grid, eta_star_val = eta_star)
	}
	
	if(process == "FARMA")
	{
		eta_star_MA = funMAMat(refinement = no_grid, eta_star_val = eta_star)	
		sim_X = funARMat(refinement = no_grid, eta_star_val = eta_star_MA)
    }
        
    # Step 4 (last n observations)
    
    sim_X_record = sim_X[,(ncol(sim_X)-(n-1)):ncol(sim_X)]
    return(sim_X_record)
}    

