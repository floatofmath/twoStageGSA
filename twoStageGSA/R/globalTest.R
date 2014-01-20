## Global test of Goeman et al 

globalTest <- function(set,perm,d,l){
  out <- gt(l,t(d[set,]),model='logistic',permutations=perm)
  p.value(out)
}
