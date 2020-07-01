# defines rank one operator defined as: (g1 \otimes g2)(f) = <f,g2>g1
rank_one_tensor <- function(g1, g2, f){
  grid <- seq(0,1,length.out=length(f))
  return( mean( f * g2(grid) ) * g1(grid) )
}
