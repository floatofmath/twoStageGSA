## Apply focus step of procedure

focus <- function(data,labels,genesets,sigSets,B=100,test="wilcoxon",side='abs',adjust="WY",reuse.p=F){
  focusWY <- function(idx){
    idx <- which(rownames(data) %in% genesets[[idx]])
    out <- mt.minP(data[idx,],labels,test = test,  B=B, side=side)
    ans <- out$adjp
    names(out) <- out$index
  }
  focusBonfH <- function(idx){
    testfun <- match.fun(test)
    idx <- rownames(data)[which(rownames(data) %in% genesets[[idx]])]
    if(reuse.p){
      return(testfun(idx))
    }
    out <- apply(data[idx,],1,function(x) testfun(x~labels)$p.value)
    ans <- p.adjust(out,method='holm')
    return(ans)
  }
  focusTest <- function(idx) {
    switch(adjust,
           WY = focusWY(idx),
           holm = focusBonfH(idx))
  }
  pWY <- lapply(sigSets,focusTest)
  return(pWY)
}
