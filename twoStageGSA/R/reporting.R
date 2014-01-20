reportScreening <- function(scr,pw.annotation,pw.kv){
  if(length(scr$sigSets)<1) return(NULL)
  if(length(scr$sigSets)==1) {
    tout <- c(scr$sigSets,select(pw.annotation,keys=scr$sigSets,keytype=pw.kv[1],columns=pw.kv[2])[,2],round(scr$qVal[scr$sigSets],3))
    names(tout) <- c('Gene Set ID','Gene Set Description','FDR adj. q-Value')
    return(tout)
  }
  tout <- data.frame(k=scr$sigSets,v=select(pw.annotation,keys=scr$sigSets,keytype=pw.kv[1],columns=pw.kv[2])[,2],e=round(scr$qVal[scr$sigSets],3))
  colnames(tout) <- c('Gene Set ID','Gene Set Description','FDR adj. q-Value')
  o <- order(tout[,3])
  rownames(tout) <- 1:nrow(tout)
  tout[o,]
}



reportFocus <- function(scr,pval,fc,q,linebreak=75,pw.annotation,gene.annotation,pw.kv,gene.kv){
  parser <- function(s){
    if(!all(s>=q)){
      idx <- which(s<q)
      ans <- paste(select(gene.annotation,keys=names(s)[idx],columns=gene.kv[2],keytype=gene.kv[1])[,2],' (',round(s[idx],3),', ',c('down','NA','up')[2+sign(fc[names(s)[idx]])],')',sep='',collapse='; ')
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
  desc <- select(pw.annotation,keys=names(out),columns=pw.kv[2],keytype=pw.kv[1])[,2]
  outtab <- cbind(names(out),desc,unlist(out))
  rownames(outtab) <- 1:length(out)
  colnames(outtab) <- c('Gene Set ID','Gene Set Description',paste('Considerable genes with local p-value <',q))
  return(outtab)
}
