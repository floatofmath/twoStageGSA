## test screening
library(twoStageGSA)
generate.dataset <- function(chips,genes,sets,setsize,nSigGenes,filt=1,rLabels=FALSE,rSets=FALSE,rSigGenes=FALSE){
    if(rLabels){
        ## random group sizes
        control <- sample(3:(chips-2),1)
        labels <- rep(0:1,c(control,chips-control))
    } else {
        if((chips %% 2 == 1)){
            stop("Uneven sample size, if rLabels = FALSE sample size must be even")
        }
        ## balanced sample sizes
        labels <- rep(0:1,each=chips/2)
    }
    ## generate genenames
    allgenes <- paste('g',10^(floor(log(genes,base=10)+1))+1:(genes*filt),sep='_')
    ## suppose that the genenames on the chip are only a subset of those
    ## that may appear in a geneset
    genenames <- sort(sample(allgenes,genes))

    if(rSets){
        ## we draw from sets from allgenes
        setlist <- lapply(1:sets,function(i){
            size <- sample(2:floor(genes*filt/4),1)
            sample(allgenes,size)
        })
    } else {
        size <- floor(length(allgenes)/sets)
        setlist <- lapply(0:(sets-1),function(i) paste('g',10^(floor(log(genes,base=10)+1))+1:(genes*filt),sep='_'))
    }
    names(setlist) <- paste('s',1:sets,sep='_')

    if(rSigGenes){
        nSigGenes <- sample(1:genes,1)
    } 
    sigGenes <- sample(genes,nSigGenes)
    data <- matrix(rnorm(genes*chips),ncol=chips)
    data[sigGenes,] <- t(t(data[sigGenes,])+labels*1000)
    return(list('data'=data,'labels'=labels,'setlist'=setlist,'genenames'=genenames,'sigGenes'=sigGenes))
}

test.screening <- function(chips,genes,sets,setsize,nSigGenes,
                           filt=1,
                           rLabels=FALSE,
                           rSets=FALSE,
                           rSigGenes=FALSE,
                           B=0,q=.05,settest='globalTest',min=6,max=100){
    simexp <- generate.dataset(chips,genes,sets,setsize,nSigGenes,filt,rLabels,rSets,rSigGenes)


    
    toSmall <- which(sapply(simexp$setlist,function(set) sum(simexp$genenames %in% set)<=min))
    toLarge <- which(sapply(simexp$setlist,function(set) sum(simexp$genenames %in% set)>=max))


    if (length(c(toSmall,toLarge))==sets){
        sigSets <- list()
    } else if(length(c(toSmall,toLarge)>0)){
        sigSets <- names(which(sapply(simexp$setlist[-c(toSmall,toLarge)],function(set) any(simexp$genenames[simexp$sigGenes] %in% set))))
    } else {
       sigSets <- names(which(sapply(simexp$setlist,function(set) any(simexp$genenames[simexp$sigGenes] %in% set))))
    }
    
    scr <- screening(simexp$data,simexp$labels==1,simexp$setlist,simexp$genenames,B=0,q=.05,settest='globalTest',min=6,max=100)
    logger <- numeric()
    if(length(sigSets) == 0){
        out <- (length(scr$sigSets) == 0)
    } else if(length(scr$sigSets) == length(sigSets)){
        out <- all(scr$sigSets %in% sigSets)
    } else {
        out <- FALSE
    }
    if(!out){
        logger <- scr$qVal[!(scr$sigSets %in% sigSets)]
    }
    ##:ess-bp-start::logger@logger:##
.ess_log_eval('logger')##:ess-bp-end:##
    
    return(out)
}

test.screening(10,100,100,30,30,filt=5,T,T,T)
sum(replicate(10000,test.screening(10,1000,20,10,20,filt=5,T,T,T)))

## two sig sets



## > 2 sig sets
