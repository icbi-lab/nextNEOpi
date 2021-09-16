#!/usr/bin/env python

"""
Check BAM file(s) if they come from PE or SE sequencing libraries

Requirements:
    * Python >= 3.6.2
    * pysam >= 0.16.0.1

Copyright (c) 202^ Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

import pysam
import sys


def check_pe(bam):
    pe = False
    read_nr = 0
    max_reads = 10000

    samfile = pysam.AlignmentFile(bam, mode="rb", threads=2)

    for read in samfile.fetch(until_eof=True):
        if read.is_paired:
            pe = True
            break
        if read_nr > max_reads:
            break
        read_nr += 1


    samfile.close()

    return pe


pe = 0
se = 0
nf = 0

for bam in sys.argv[1:]:
    pe_tmp = check_pe(bam)
    if pe_tmp:
       pe += 1
    else:
       se += 1
    nf += 1

if se == nf:
    print("SE", end="")
elif pe == nf:
    print("PE", end="")
else:
    print("MIXED", end="")
