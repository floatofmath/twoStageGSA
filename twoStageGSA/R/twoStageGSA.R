## Private functions

fdr <- function(pv,al=.05){
  ## maximum p-value below fdr alpha
  if (sum(is.na(pv))>0 | sum(pv<0,na.rm=T)>0 | sum(pv>1,na.rm=T)>0) {
    warning("Some p-values are not-valid or missing! These p-values are ignored")
  }
  spv = sort(pv[!is.na(pv)])
  G = length(spv)
  temp = spv<=seq(1,G,1)*al/G
  pID = ifelse(length(spv[temp==T])>0,max(spv[temp==T]),NA)
  pID
} 

