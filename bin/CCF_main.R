source("CCF_functions.R")

exdata<-read.csv("ex_data.tsv", sep="\t",
  stringsAsFactors=TRUE)

CCFest<-getCCF(exdata)

