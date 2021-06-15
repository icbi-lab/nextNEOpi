#!/bin/sh

# gencode.v33.primary_assembly.annotation.gtf
GTF=$1
cut -f1-5 $1 | grep -w "exon" | grep -E "^chr" | cut -f1,4,5 | bedtools sort | bedtools merge  > gencode.v33.primary_assembly.annotation.exon_merged.bed