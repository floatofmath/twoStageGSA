## Global test based on inverse normal combination procedure

screening_invnormT <- function(set,perm,d,l){
  p <- apply(d[set,],1,function(y) t.test(y~l)$p.value)
  q <- sum(qnorm(1-p))/sqrt(length(p))
  return(1-pnorm(q))
}
