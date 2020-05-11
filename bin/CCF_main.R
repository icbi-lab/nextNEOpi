source("CCF_functions.R")

args = commandArgs(trailingOnly=TRUE)

varDataIn <- args[1]
outFile <- args[2]

varData<-read.csv(varDataIn, sep="\t", stringsAsFactors=TRUE)

CCFest<-getCCF(varData)

write.table(CCFest, file=outFile, sep="\t", quote=FALSE)
