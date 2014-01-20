## Hotelling's t-test

hotellingT <- function(set,perm,d,l){
  hotelling <- function(d1,d2){
    k <- ncol(d1)
    n1 <- nrow(d1)
    n2 <- nrow(d2)
    if(k > n1 | k > n2) {
      mysolve <- ginv
    } else {
      mysolve <- solve
    }
    xbar1 <- apply(d1,2,mean)
    xbar2 <- apply(d2,2,mean)
    dbar <- xbar2-xbar1
    v <- ((n1-1)*var(d1)+(n2-1)*var(d2))/(n1+n2-2)
    t2 <- n1*n2*dbar%*%mysolve(v)%*%dbar/(n1+n2)
    f <- (n1+n2-k-1)*t2/((n1+n2-2)*k)
    return(1-pf(f,k,n1+n2-k-1))
  }
  if(is.factor(l)){
    ls <- levels(l)
  } else {
    ls <- unique(l)
  }
  if(length(ls)>2) stop('Too many groups!')
  d1 <- d[set,l==ls[1]]
  d2 <- d[set,l==ls[2]]
  hotelling(t(d1),t(d2))
}
