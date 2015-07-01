## Apply focus step of procedure







#' Perform focus step of two-stage GSA
#' 
#' Performs focus step of two-stage GSA. In this step single genes are tested
#' at significance levels resulting from the screening step, which only tests
#' gene sets.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param data Matrix containing expression values. Genes are in rows, samples
#' in columns.
#' @param labels Vector containing gene identifiers that correspond to those in
#' genesets
#' @param genesets List of genesets, given as vectors of gene identifiers
#' @param sigSets Output from screening step
#' @param genes Vector with gene identifier
#' @param B Number of resampling iterations (e.g. for Westfall & Young)
#' @param test name of test function, if reuse.p is FALSE a test function that
#' given if reuse.p is TRUE a function that given takes a formula argument
#' (x~labels) and returns an object with $p.value slot. Vector of gene
#' identifiers returns the corresponding unadjusted p-values
#' @param side should one or two-sided tests be performed
#' @param adjust multiplicity adjustment to use
#' @param reuse.p Sometimes it is desired to re-use the p-values of a previous
#' differential gene expression analysis. If TRUE test should be a vector of
#' unadjusted p-values
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (data, labels, genesets, sigSets, B = 100, test = "wilcoxon", 
#'     side = "abs", adjust = "WY", reuse.p = F) 
#' {
#'     focusWY <- function(idx) {
#'         idx <- which(rownames(data) %in% genesets[[idx]])
#'         out <- mt.minP(data[idx, ], labels, test = test, B = B, 
#'             side = side)
#'         ans <- out$adjp
#'         names(out) <- out$index
#'     }
#'     focusBonfH <- function(idx) {
#'         testfun <- match.fun(test)
#'         idx <- rownames(data)[which(rownames(data) %in% genesets[[idx]])]
#'         if (reuse.p) {
#'             return(testfun(idx))
#'         }
#'         out <- apply(data[idx, ], 1, testfun)
#'         ans <- p.adjust(out, method = "holm")
#'         return(ans)
#'     }
#'     focusTest <- function(idx) {
#'         switch(adjust, WY = focusWY(idx), holm = focusBonfH(idx))
#'     }
#'     pWY <- lapply(sigSets, focusTest)
#'     return(pWY)
#'   }
#' 
#' @export focus
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
