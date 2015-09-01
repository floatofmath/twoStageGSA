##' Various global test procedures that can be applied in the screening stage of twoStageGSA
##'
##' @title Global test procedures 
##' @param set identifiers of probes included in the gene set
##' @param perm number of permutations to be applied in resampling tests
##' @param d numeric matrix containing expression values
##' @param l vector or factor containing treatment/phenotype labels
##' @return p-value of the test applied to the data 
##' @author Florian Klinglmueller
##' @name screening-tests
NULL

##' @rdname screening-tests
screening_hotellingT <- function(set,perm,d,l){
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

##' @rdname screening-tests
screening_globalTest <- function(set,perm,d,l){
  if(is.null(rownames(d))) {
      ## somehow gt does require rownames
      rownames(d) <- 1:nrow(d)
  }
  out <- gt(l,t(d[set,]),model='logistic',permutations=perm)
  p <- p.value(out)
  return(p)
}

##' @rdname screening-tests
screening_invnormT <- function(set,perm,d,l){
  p <- apply(d[set,],1,function(y) t.test(y~l)$p.value)
  q <- sum(qnorm(1-p))/sqrt(length(p))
  return(1-pnorm(q))
}

##' @rdname screening-tests
screening_nettleton <- function(set,perm,d,l) {
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

  mrpp(l,e.com(t(d[set,])),nperm=perm)
}

