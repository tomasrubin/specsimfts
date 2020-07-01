# defines a compact operator through its kernel "ker" applied on the function "f"
kernel_operator <- function( ker, f ){
  grid <- seq(0,1,length.out=length(f))
  grid_matrix <- kronecker(grid,matrix(1,1,n_grid))
  ker( grid_matrix, t(grid_matrix) ) %*% f / length(grid)
}