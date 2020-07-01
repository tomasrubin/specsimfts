###################################################################################################
## create kernel of an integral operator given its eigendecomposition
kernel_from_eig <- function( eigenvalues, eigenfunctions ){
  return(function(x,y){
    r <- 0
    for (ii in 1:length(eigenvalues)){
      r <- r + eigenvalues[ii] * eigenfunctions[[ii]](x) * eigenfunctions[[ii]](y)
    }
    return(r)
  })
}