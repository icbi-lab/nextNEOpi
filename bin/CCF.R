args = commandArgs(trailingOnly=TRUE)

varDataIn <- args[1]
outFile <- args[2]

varData<-read.csv(varDataIn, sep="\t", stringsAsFactors=TRUE)


getCCF<-function(data) {
  
  # Restrict to Chromosomes 1-22
  data<-data[which(as.character(data$Chr) %in% 
      as.character(seq(1:22))),]
  
  # Select wanted columns
  wantcol<-c("SubjectID",
    "SampleType", 
    "Chr",
    "pos",
    "normCnt",
    "mutCnt",
    "VAF",
    "ploidy",
    "cnA",
    "cnB",
    "purity")
  ind<-match(wantcol, colnames(data))
  if (sum(is.na(ind)>0)) stop("Wrong data format!\n")
  indata<-data[,ind]
  
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
  local.copy.number<-idata[,"cnA"]
  
  x<-dbinom(n.alt, depth, 
    prob=sapply(seq(0.01,1,length.out=100),
      f.function,purity,local.copy.number))
  
  if(min(x)==0) x[length(x)]<-1
  
  names(x)<-seq(0.01,1,length.out=100)
  est<-sub.cint(x, n.alt, depth)
  
  return(est)
  
}

CCFest<-getCCF(varData)

write.table(CCFest, file=outFile, sep="\t", quote=FALSE)
