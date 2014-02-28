## Global test of Goeman et al 

screening_globalTest <- function(set,perm,d,l){
  if(is.null(rownames(d))) {
      ## somehow gt does require rownames
      rownames(d) <- 1:nrow(d)
  }
  out <- gt(l,t(d[set,]),model='logistic',permutations=perm)
  p <- p.value(out)
  return(p)
}

