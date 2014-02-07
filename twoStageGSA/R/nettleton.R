## Implements Nettletons (resampling based) gene set test
################################################################################
# Author of the following 3 functions:  Dan Nettleton                          #
# Date: August 7, 2007                                                         #
# Downloaded from: http://www.public.iastate.edu/~dnett/useR/GeneCatTesting.txt#
################################################################################
screening_nettleton <- function(set,perm,d,l) {
  mrpp(l,e.com(t(d[set,])),nperm=perm)
}

e.com=function (y) {#Nettleton
  K=ncol(y)
  scale=rep(0,K)
  for(i in 1:K){
    D=dist(y[,i])
    scale[i]=1/sum(D)
  }
  y=sweep(y,2,scale,FUN="*")
  return(y)
}

get.delta=function (x,D,n,N,w) {#Nettleton
  inter.sum=tapply(1:N,x,FUN="inter.dist",D)
  inter.mean=inter.sum*2/(n*(n-1))
  delta=sum(w*inter.mean)
  return(delta)
}

inter.dist=function (index,D) {#Nettleton
  return(sum(as.matrix(D)[index,index])/2)
}

mrpp=function (x,y,nperm) { #Nettleton's p-value for a gene set. 
  D=dist(y)
  x=as.factor(x)
  n=table(x)
  N=sum(n)
  w=n/N
  delta.obs=get.delta(x,D,n,N,w)
  delta <- unlist(mylapply(1:nperm,function(i) {get.delta(sample(x),D,n,N,w)}))
  return(p=mean(delta.obs>=delta))
}
