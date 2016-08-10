## Apply screening stage of procedure

#' Performs the screening step
#' 
#' twoStageGSA provides a couple of screening test procedures. Those
#' are 'globalTest', 'invNormal', 'hotellingT', and
#' 'nettleton'. Custom functions may be used. In that case settest has
#' to point to a function with the name matching
#' paste('screening_',settest). The corresponding function needs to
#' return an unadjusted p-value for a given geneset. See
#' \code{\link{screening-tests}} and the vignette for further details
#' and some examples.
#' 
#' @param data Expression values as a matrix with genes in rows
#' @param labels Object defining the group labels of the dataset. E.g. for
#' globalTest this may be a factor.
#' @param genesets A list whose elements are vectors of gene identifiers
#' representing the gene sets.
#' @param genes vector of gene identifiers specifying the gene each row of data
#' corresponds to and that match to the identifyers used to define gene sets.
#' If not given rownames are assumed to be gene identifyers.
#' @param B Number of randomizations to use if a randomization test is used
#' @param q q-Value cutoff for selection
#' @param settest A character string giving the name of the screening test
#' procedure. Either one of the methods already implemented in the package, or
#' a custom test function (see Details).
#' 
#' @param min minimum number of genes per set for which measurements are available
#' @param max maximum number of genes per set for which measurements are available
#' @param parcomp should parallel processing be used and which implementation of mclapply 
#' @export screening
screening <- function(data,labels,genesets,genes=rownames(data),B=0,q=.05,settest='globalTest',min=6,max=100,parcomp=c('none','parallel','bt88.03.704')){
  if(is.null(names(genesets))){
      stop("names(genesets) is NULL: genesets must have names!")
  }
  ## prune genes which are not part of any geneset
  toSmall <- which(sapply(genesets,function(set) sum(genes %in% set)<=min))
  toLarge <- which(sapply(genesets,function(set) sum(genes %in% set)>=max))
  if(length(c(toSmall,toLarge))>0){
      genesets <- genesets[-c(toSmall,toLarge)]
  }
  if(!length(genesets)>=1){
      warning("No sets that have more than minimum required/less than maximum allowed members in dataset")
      return(list(sigSets=character(),qVal=numeric()))
  }
  # remove rows that are in no geneset
  inSet <- unlist(sapply(genesets,function(set) which(genes %in% set)))
  #impW <- tabulate(inSet)/length(inSet)
  inSet <- sort(unique(inSet))
  data <- data[inSet,]
  genes <- genes[inSet]
                                        # in case not all rows are unique genes
  idSets <- lapply(genesets,function(set) which(genes %in% set))
  ## some helpful constants
  S <- length(genesets)
  L <- sapply(genesets,length)
  n <- max(L)
  G <- nrow(data)
  m <- ncol(data)
  
  ## compute set test p.value and correct FDR
  parcomp <- match.arg(parcomp)
  myapply <- switch(parcomp,
                    'none'=lapply,
                    'parallel'=mclapply,
                    'bt88.03.704'=bt88.03.704::mclapply2)
  fun <- match.fun(paste("screening_",settest,sep=''))
  pmrpp <- simplify2array(myapply(idSets,function(set) fun(set,perm=B,d=data,l=labels)))
  names(pmrpp) <- names(genesets)
  q.vals <- p.adjust(pmrpp,method='fdr')

  ##################### sets significant at the alpha FDR level
  sigSets <- names(genesets)[which(q.vals <= q)]
  return(list(sigSets=sigSets,qVal=q.vals,size=L))

}

