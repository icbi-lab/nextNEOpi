#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "./tpm.out.txt"
  args[3] = "./rpkm.out.txt"
}

# Import count matrix
countdata <- read.table(args[1], header = TRUE)
counts <- data.frame(countdata[1], countdata[6], countdata[7])
names(counts) <- c("Geneid", "Length", "Counts")

# Define TPM/RPKM functions
mytpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

myrpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

# Calulate TPM/RPKM and write output
count_tpm <- mytpm(counts$Counts, counts$Length)
count_rpkm <- myrpkm(counts$Counts, counts$Length)

out1 <- data.frame(counts$Geneid, counts$Counts, count_tpm)
names(out1) <- c('Geneid', 'Raw_Counts', 'TPM')
write.table(out1, file = args[2], sep = "\t", row.names = FALSE, quote = FALSE)

out2 <- data.frame(counts$Geneid, counts$Counts, count_rpkm)
names(out2) <- c('Geneid', 'Raw_Counts', 'RPKM')
write.table(out2, file = args[3], sep = "\t", row.names = FALSE, quote = FALSE)
