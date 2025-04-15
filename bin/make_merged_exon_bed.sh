#!/bin/sh

# gencode.v47.primary_assembly.annotation.gtf
GTF=$1
cut -f1-5 $GTF | grep -w "exon" | grep -E "^chr" | cut -f1,4,5 | bedtools sort | bedtools merge  > gencode.v47.primary_assembly.annotation.exon_merged.bed
