## Apply focus step of procedure

focus <- function(data,labels,genesets,sigSets,genes=rownames(data),B=100,test="wilcoxon",side='abs',adjust="WY",reuse.p=F){
  if(reuse.p){
      if(!(length(genes) == length(test))){
          stop(paste("Number of genes (",length(genes),") needs to match number of p-values (",length(test),")",sep=''))
      }
  }
  focusWY <- function(idx){
      if(reuse.p){
          stop("Re-using p-values with Westfall & Young is not possible")
      }
    gidx <- which(genes %in% genesets[[idx]])
    out <- mt.minP(data[gidx,],labels,test = test,  B=B, side=side)
    ans <- out$adjp
    names(out) <- out$index
  }
  focusBonfH <- function(idx){
      if(reuse.p){
          
          gidx <- which(genes %in% genesets[[idx]])
          
          out <- test[gidx]
      } else {
          testfun <- match.fun(test)
          gidx <- which(genes %in% genesets[[idx]])
          out <- apply(data[gidx,],1,function(x) testfun(x~labels)$p.value)
      }
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
