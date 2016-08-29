#' Produces human readable summary report
#' 
#' @param scr Output object from screening step of twoStageGSA
#' @param pw.annotation Annotation environment that maps pathway identifiers to
#' some annotation
#' @param pw.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' @author Florian Klinglmueller 
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
#' 
#' @param scr output from screening step
#' @param pval 
#' @param fc 
#' @param q 
#' @param linebreak 
#' @param pw.annotation Annotation environment that maps pathway identifiers to
#' some annotation
#' @param gene.annotation Annotation environment that maps gene identifiers to
#' some annotation
#' @param pw.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' @param gene.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' w
#' @author Florian Klinglmueller
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
      ## names(pval) <- sets[[scr$sigSets]]
      return(parser(unlist(pval)))
  }
  ## names(pval) <- scr$sigSets
  ## pval <- lapply(names(pval),function(set) {
  ##                    ans <- pval[[set]]
  ##                    names(ans) <- sets[[set]]
  ##                    ans
  ##                })
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


#' Produces a nice plot 
#' 
#' @param rank rank of set to plot
#' @param scr output from screening step
#' @param focus focus level result
#' @param fc vector of fold-changes/test statistics
#' @param q local significance level
#' @param linebreak 
#' @param pw.annotation Annotation environment that maps pathway identifiers to
#' some annotation
#' @param gene.annotation Annotation environment that maps gene identifiers to
#' some annotation
#' @param pw.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' @param gene.kv Character value of the form c('KEY','VALUE') specifying which
#' key- and valuetypes should be used in the corresponding annotation database
#' w
#' @author Florian Klinglmueller
#' 
#' @export reportFocus
plot_2SGSAS <- function(scr,focus,fc,q,linebreak=75,pw.annotation,gene.annotation,pw.kv,gene.kv){
  parser <- function(s){
    if(!all(s>=q)){
      idx <- which(s<q)
      
      anno <- AnnotationDbi:::select(gene.annotation,keys=names(s),columns=gene.kv[2],keytype=gene.kv[1])      
      annotation <- join.duplicated(anno)
      ans <- paste(annotation,' (',round(s,3),', ',c('down','NA','up')[2+sign(fc[names(s)])],')',sep='',collapse='; ')      
      
      
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
    return(parser(unlist(focus)))
  }
  
  out <- sapply(focus,parser)
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
