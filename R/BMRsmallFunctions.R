fun_hc <- function(n1,N1,n2,N2)
{
    p = 0
    for ( ni in 0:n1){
        p = p+ fun_hp(ni,N1,n2,N2)
    }
    p[p > 1] <- 1
    p[p < 0] <- 0
    return(p)
}

fun_hp <- function(n1,N1,n2,N2)
{
  p <- exp(lgamma(N1+1)+lgamma(N2+2)+lgamma(n1+n2+1)+lgamma(N1+N2-n1-n2+1)-
             lgamma(n1+1)-lgamma(n2+1)-lgamma(N1-n1+1)-lgamma(N2-n2+1)-
             lgamma(N1+N2+2))
  p[p > 1] <- 1
  p[p < 0] <- 0
  return(p)
}
