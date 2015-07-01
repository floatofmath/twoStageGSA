#' Produces human readable summary report
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param scr %% ~~Describe \code{scr} here~~
#' @param pw.annotation Annotation environment that maps pathway identifiers to
#' some annotation
#' @param pw.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
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
#' function (scr) 
#' {
#'     if (length(scr$sigSets) < 1) 
#'         return(NULL)
#'     if (length(scr$sigSets) == 1) {
#'         tout <- c(scr$sigSets, round(scr$qVal[scr$sigSets], 3), 
#'             unlist(mget(sub("hsa", "", scr$sigSets), KEGGPATHID2NAME, 
#'                 ifnotfound = NA)))
#'         names(tout) <- c("Gene Set ID", "Gene Set Description", 
#'             "FDR adj. q-Value")
#'         return(tout)
#'     }
#'     tout <- cbind(scr$sigSets, unlist(mget(sub("hsa", "", scr$sigSets), 
#'         KEGGPATHID2NAME, ifnotfound = NA)), round(scr$qVal[scr$sigSets], 
#'         3))
#'     tout <- tout[order(tout[, 2]), ]
#'     colnames(tout) <- c("Gene Set ID", "Gene Set Description", 
#'         "FDR adj. q-Value")
#'     o <- order(tout[, 3])
#'     rownames(tout) <- 1:nrow(tout)
#'     tout[o, ]
#'   }
#' 
#' @export reportScreening
reportScreening <- function(scr,pw.annotation,pw.kv){
  if(length(scr$sigSets)<1) return(NULL)
  annotation <- join.duplicated(AnnotationDbi:::select(pw.annotation,keys=scr$sigSets,keytype=pw.kv[1],columns=pw.kv[2]))
  if(length(scr$sigSets)==1) {
    tout <- c(scr$sigSets,annotation,scr$size[scr$sigSets],round(scr$qVal[scr$sigSets],3))
    names(tout) <- c('Gene Set ID','Gene Set Description','# Genes','FDR adj. q-Value')
    return(tout)
  }
  tout <- data.frame(k=scr$sigSets,v=annotation,l=scr$size[scr$sigSets],e=round(scr$qVal[scr$sigSets],3))
  colnames(tout) <- c('Gene Set ID','Gene Set Description','# Genes','FDR adj. q-Value')
  o <- order(tout[,4])
  rownames(tout) <- 1:nrow(tout)
  tout[o,]
}










#' Produces a nice report
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param scr %% ~~Describe \code{scr} here~~
#' @param pval %% ~~Describe \code{pval} here~~
#' @param fc %% ~~Describe \code{fc} here~~
#' @param q %% ~~Describe \code{q} here~~
#' @param linebreak %% ~~Describe \code{linebreak} here~~
#' @param pw.annotation Annotation environment that maps pathway identifiers to
#' some annotation
#' @param gene.annotation Annotation environment that maps gene identifiers to
#' some annotation
#' @param pw.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' @param gene.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' w
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
#' function (scr, pval, fc, q, linebreak = 75) 
#' {
#'     parser <- function(s) {
#'         if (!all(s >= q)) {
#'             idx <- which(s < q)
#'             ans <- paste(mget(paste(names(s)[idx], "at", sep = "_"), 
#'                 hugene10stv1hsentrezgSYMBOL), " (", round(s[idx], 
#'                 3), ", ", c("down", "NA", "up")[2 + sign(fc[paste(names(s)[idx], 
#'                 "at", sep = "_"), ])], ")", sep = "", collapse = "; ")
#'             if (nchar(ans) > linebreak) {
#'                 splt <- c(seq(1, nchar(ans), linebreak), nchar(ans))
#'                 ans <- paste(strwrap(ans, linebreak), collapse = "\n")
#'             }
#'         }
#'         else {
#'             return(NULL)
#'         }
#'     }
#'     if (length(scr$sigSets) < 1) 
#'         return(NULL)
#'     if (length(scr$sigSets) == 1) {
#'         return(parser(unlist(pval)))
#'     }
#'     out <- sapply(pval, parser)
#'     names(out) <- scr$sigSets
#'     out <- out[!sapply(out, is.null)]
#'     desc <- unlist(mget(sub("hsa", "", names(out)), KEGGPATHID2NAME))
#'     outtab <- cbind(names(out), desc, unlist(out))
#'     rownames(outtab) <- 1:length(out)
#'     colnames(outtab) <- c("Gene Set ID", "Gene Set Description", 
#'         paste("Considerable genes with local p-value <", q))
#'     return(outtab)
#'   }
#' 
#' @export reportFocus
reportFocus <- function(scr,pval,fc,q,linebreak=75,pw.annotation,gene.annotation,pw.kv,gene.kv){
  parser <- function(s){
    if(!all(s>=q)){
      idx <- which(s<q)
      
      anno <- AnnotationDbi:::select(gene.annotation,keys=names(s)[idx],columns=gene.kv[2],keytype=gene.kv[1])      
      annotation <- join.duplicated(anno)
      ans <- paste(annotation,' (',round(s[idx],3),', ',c('down','NA','up')[2+sign(fc[names(s)[idx]])],')',sep='',collapse='; ')      
      
      
      if(nchar(ans)>linebreak){
        splt <- c(seq(1,nchar(ans),linebreak),nchar(ans))
        ans <- paste(strwrap(ans,linebreak),collapse='\n')
      }
      return(ans)
    } else {
      return(NULL)
    }
  }
  if(length(scr$sigSets)<1) return(NULL)
  if(length(scr$sigSets)==1){
    return(parser(unlist(pval)))
  }
  
  out <- sapply(pval,parser)
  names(out) <- scr$sigSets
  out <- out[!sapply(out,is.null)]
  desc <- join.duplicated(AnnotationDbi:::select(pw.annotation,keys=names(out),columns=pw.kv[2],keytype=pw.kv[1]))
  
  outtab <- cbind(names(out),desc,unlist(out))
  rownames(outtab) <- 1:length(out)
  colnames(outtab) <- c('Gene Set ID','Gene Set Description',paste('Considerable genes with local p-value <',q))
  return(outtab)
}

join.duplicated <- function(annotation,return.column=2){
  if(anyDuplicated(annotation[,1])){
      annotation <- tapply(annotation[,return.column],annotation[,1],paste,sep=' ',collapse='|')
  } else {
      annotation <- annotation[,return.column]
  }
  return(annotation)
}

## join.duplicated <- function(annotation,key=1){
##     joiner <- function(char){
##         paste(char,sep='|')
##     }
##     ukeys <- annotation[,key]
##     if(!anyDuplicated(ukeys)){
##         return(annotation)
##     } else if(ncol(annotation) > 2){
##         anno <- apply(annotation[,-key],2,function(c) unlist(tapply(c,ukeys,joiner)))
##     } else if(ncol(annotation) == 2){
##         anno <- unlist(tapply(annotation[,-key],ukeys,joiner))
##     }
##     out <- cbind(rownames(anno),anno)
##     colnames(out) <- colnames(annotation)
##     return(out)    
## }
