##' Searches for a string in the title and description of a number of pathway databases
##'
##' @title Search related pathways
##' @param string character string to search for
##' @return 
##' @author float
search.pw <- function(string){
  n <- grep(string,nci,ignore.case=T)
  r <- grep(string,reactome,ignore.case=T)
  k <- grep(string,kegg.names,ignore.case=T)
  b <- grep(string,broad.names,ignore.case=T)
  out <- list()
  if(length(n)>0){
    names(n)<- nci[n]
    out$nci <- n
    out$nci.length <- sapply(NCI.cyList[n],numNodes)
  } 
  if(length(r)>0){
    names(r) <- reactome[r]
    out$reactome <- r
  }
  if(length(k)>0){
    names(k) <- kegg.names[k]
    out$kegg <- k
  }
  if(length(b)>0){
    names(b) <- broad.names[b]
    out$broad <- b
  }
  return(out)
}

##' Extracts uniprot IDS from NCI data
##'
##' @title Get uniprot ID
##' @param ID pathway ID
##' @param cyList 
##' @return 
##' @author float
getUniprot <- function(ID,cyList=NCI.cyList){
  try(nd <- unlist(nodeData(NCI.cyList[[ID]],attr='biopax.xref.UNIPROT')),silent=T)
  if(exists('nd')) unique(nd[nd!='unassigned'])
  else NULL
}
##' Extracts entrezgene IDS from NCI data
##'
##' @title Get entrezgeneID
##' @inheritParams getUniprot
##' @return 
##' @author float
getEntrezgene <- function(ID,cyList=NCI.cyList){
  try(nd <- unlist(nodeData(NCI.cyList[[ID]],attr='biopax.xref.ENTREZGENE')),silent=T)
  if(exists('nd')) unique(nd[nd!='unassigned'])
  else NULL
}

##' Extract affy probe ids from cyList pathways
##'
##' @title cyList to affy IDs
##' @param ID pathway ID
##' @param db pathway database
##' @param ... further arguments passed to getUniprot and getEntrezgene
##' @return 
##' @author float
cylist2afid <- function(ID,env='hugene10stv1hsentrezg',...){
  up <- getUniprot(ID,...)
  eg <- getEntrezgene(ID,...)
  upmap <- eval(substitute(revmap(foo), list(foo=as.name(paste(env,'UNIPROT',sep='')))))
  egmap <- eval(substitute(revmap(foo), list(foo=as.name(paste(env,'ENTREZID',sep='')))))
  out <- character(0)
  if(!is.null(up))
    out <- na.omit(unlist(mget(up,envir=upmap,ifnotfound=NA)))
  if(!is.null(eg))
    out <- na.omit(c(out,unlist(mget(eg,envir=egmap,ifnotfound=NA))))
  return(out)
}

##' Extract affy probe ids from BROAD pathways
##'
##' @title BROAD to affy IDs
##' @param ID pathway ID
##' @param db pathway database
##' @param env annotation environment that maps Gene Symbols to probids
##' @return 
##' @author float
broad2afid <- function(ID,db=broad,env=hugene10stv1hsentrezgSYMBOL){
  ids <- unlist(slot(broad[[ID]],name='geneIds'))
  map <- revmap(env)
  out <- unlist(mget(ids,envir=map,ifnotfound=NA))
  attr(out,"description") <- slot(broad[[ID]],'shortDescription')
  return(out)
}

##' Matches a list of Affymetrix Probeset Ids to those available in a dataset
##'
##' @title Match Affy Ids
##' @param afids list of Affymetrix probeset IDs
##' @param data  \code{expressionSet} 
##' @param nmin number of features that should be shared at least (if less are found \code{NULL} is returned)
##' @return 
##' @author float
match.afid <- function(afids,data=data.fil,nmin){
  out <- afids[which(afids %in% featureNames(data))]
  attributes(out) <- attributes(afids[which(afids %in% featureNames(data))])
  if(length(out)<nmin){
    return(NULL)
  }
  return(out)
}
##' Turn the results of a keyword search into a pathway list
##'
##' @title Make pathway list
##' @param object Search result object
##' @param data dataset to which IDs should be matched if NULL no matching is performed
##' @param nmin minimum number of features per set
##' @return 
##' @author float
make.pws <- function(object,data=NULL,nmin=10){

  pw <- list()
  if('nci' %in% objects(object)){
    pw <- c(pw,lapply(object$nci,cylist2afid))
  }
  if('broad' %in% objects(object)){
    pw <- c(pw,lapply(object$broad,broad2afid))

  }
  if(!is.null(data)){
    pw <- lapply(pw,match.afid,nmin=nmin,data=data)
    nn <- sapply(pw,is.null)
    pw <- pw[which(!nn)]
  }
  return(pw)
}
    
##' Make a human readable list of selected pathways
##'
##' @title List pathways
##' @param object Search result object
##' @return 
##' @author float
list.pw <- function(object){
  if('nci' %in% objects(object)){
    cat("** NCI Pathways\n - ")
    cat(names(object$nci),sep='\n- ')
  }
  if('broad' %in% objects(object)){
    desc <- lapply(object$broad,function(o) attr(o,'description'))
    nd <- paste(gsub('_',' ',names(object$broad)),desc,sep=" :: ")
    cat("** BROAD Pathways\n - ")
    cat(nd,sep='\n- ')
  }
  ''
}
  
## hif1a <- search.pw('HIF')
## nfkb <- search.pw('nfk|kappa')
## nrf2 <- search.pw('nrf|erythroid')
## cmyc <- search.pw('c-myc')
## sp1 <- search.pw('sp-1')
## ap1 <- search.pw('ap-1')
## save(hif1a,nfkb,nrf2,cmyc,sp1,ap1,file='pwsearch.Rd')
