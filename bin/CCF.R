getCCF<-function(indata) {
  
  # remove chr prefix, will be added again at the end
  indata$Chr<-gsub("Chr", "", gsub("chr","", indata$Chr))
  
  # Restrict to Chromosomes 1-22
  indata<-indata[which(as.character(indata$Chr) %in% 
      as.character(c(seq(1:22), "X", "Y"))),]
  
  # Select wanted columns
  wantcol<-c("Chr",
    "pos",
    "normCnt",
    "mutCnt",
    "cnA",
    "cnB",
    "purity")
  ind<-match(wantcol, colnames(indata))
  if (sum(is.na(ind)>0)) stop("Wrong data format!\n")
  
  # Append new results
  outdata<-matrix(NA, nrow=nrow(indata), ncol=5)
  colnames(outdata)<-c("CCF",
    "CCF.05", 
    "CCF.95",
    "pSubclonal",
    "pClonal")
  
  for (i in 1:nrow(indata)) {
    
    outdata[i,]<-absCCF2(indata[i,])
    
  }
  
  outdata<-cbind(indata, outdata)
  outdata$Chr<-paste0("chr", outdata$Chr)
    
  return(outdata)
  
}

sub.cint<-function(x, n.alt, depth, prob=0.95) {
  
  xnorm<-x/sum(x)
  xsort<-sort(xnorm, decreasing=TRUE)
  xcumLik<- cumsum(xsort)
  n<-sum(xcumLik < prob)+1
  LikThresh<-xsort[n]
  cint<-x[xnorm>=LikThresh]
  all<-as.numeric(names(x))
  cellu<-as.numeric(names(cint))
  l.t<-cellu[1]
  r.t<-cellu[length(cellu)]
  m<-cellu[which.max(cint)]
  
  prob.subclonal<-sum(xnorm[1:90])
  prob.clonal<-sum(xnorm[91:100])
  
  out<-c(abs.CCF=m,
    abs.CCF.5=l.t, 
    abs.CCF.95=r.t,
    pSub=prob.subclonal,
    pClo=prob.clonal)
  
  return(out)
  
}

f.function<-function (c, purity, local.copy.number) {
 
   res<-min(c((purity*c)/(2*(1-purity)+purity*local.copy.number),1))
  
  return(res)
   
}

absCCF2<-function(idata) {
  
  n.alt<-idata[,"mutCnt"]
  depth<-idata[,"mutCnt"]+idata[,"normCnt"]
  purity<-idata[,"purity"]
  local.copy.number<-idata[,"cnA"]+idata[,"cnB"]
  
  missingValues<-sum(is.na(c(n.alt, depth, purity, local.copy.number)))>0
  if (missingValues) return(NA)
  
  x<-dbinom(n.alt, depth, 
    prob=sapply(seq(0.01,1,length.out=100),
      f.function,purity,local.copy.number))
  
  if(min(x)==0) x[length(x)]<-1
  
  names(x)<-seq(0.01,1,length.out=100)
  est<-sub.cint(x, n.alt, depth)
  
  return(est)
  
}

# Inputs
args = commandArgs(trailingOnly=TRUE)
varDataIn<-args[1]
outFile<-args[2]

# Analysis
varData<-read.table(varDataIn, sep="\t", header=TRUE)
CCFest<-getCCF(varData)

# Output
write.table(CCFest, file=outFile, sep="\t", quote=FALSE, row.names=FALSE)
