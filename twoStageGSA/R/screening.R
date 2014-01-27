## Apply screening stage of procedure

screening <- function(data,labels,genes=rownames(data),genesets,B=0,q=.05,settest='globalTest',min=6,max=100){
  ## prune genes which are not part of any geneset
  toSmall <- which(sapply(genesets,function(set) sum(genes %in% set)<=min))
  toLarge <- which(sapply(genesets,function(set) sum(genes %in% set)>=max))
  if(length(c(toSmall,toLarge))>0){
      genesets <- genesets[-c(toSmall,toLarge)]
  }
  inSet <- unlist(sapply(genesets,function(set) which(genes %in% set)))
  impW <- tabulate(inSet)/length(inSet)
  inSet <- sort(unique(inSet))
  data <- data[inSet,]
  ## some helpful constants
  S <- length(genesets)
  n <- max(sapply(genesets,length))
  G <- nrow(data)
  m <- ncol(data)

  ## compute nettleton statistic and correct FDR
  pmrpp <- sapply(genesets,match.fun(settest),perm=B,d=data,l=labels)
  names(pmrpp) <- names(genesets)
  q.vals <- p.adjust(pmrpp,method='fdr')

  ##################### sets significant at the alpha FDR level
  sigSets <- names(genesets)[which(q.vals <= q)]
  return(list(sigSets=sigSets,qVal=q.vals))

}
