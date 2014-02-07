## Apply screening stage of procedure

screening <- function(data,labels,genesets,genes=rownames(data),B=0,q=.05,settest='globalTest',min=6,max=100){
  ## prune genes which are not part of any geneset
  toSmall <- which(sapply(genesets,function(set) sum(genes %in% set)<=min))
  toLarge <- which(sapply(genesets,function(set) sum(genes %in% set)>=max))
  if(length(c(toSmall,toLarge))>0){
      genesets <- genesets[-c(toSmall,toLarge)]
  }

  # remove rows that are in no geneset
  inSet <- unlist(sapply(genesets,function(set) which(genes %in% set)))
  impW <- tabulate(inSet)/length(inSet)
  inSet <- sort(unique(inSet))
  data <- data[inSet,]

  # in case not all rows are unique genes
  IDsets <- lapply(genesets,function(set) which(genes %in% set))
  ## some helpful constants
  S <- length(genesets)
  n <- max(sapply(genesets,length))
  G <- nrow(data)
  m <- ncol(data)

  ## compute set test p.value and correct FDR
  pmrpp <- sapply(IDsets,match.fun(paste("screening_",settest,sep='')),perm=B,d=data,l=labels)
  names(pmrpp) <- names(genesets)
  q.vals <- p.adjust(pmrpp,method='fdr')

  ##################### sets significant at the alpha FDR level
  sigSets <- names(genesets)[which(q.vals <= q)]
  return(list(sigSets=sigSets,qVal=q.vals))

}
