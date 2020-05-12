#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2
    * Pandas
    * NumPy

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '3', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''

import argparse
# import math
import os,sys
import pandas as pd
# import numpy as np

def append_score(pvacseq_tsv, immun_tsv, output):
    pvacseq_reader = pd.read_csv(pvacseq_tsv, sep = '\t')
    pvacseq_df = pd.DataFrame(pvacseq_reader)
    print(pvacseq_df)
    imm_reader = pd.read_csv(immun_tsv, sep = '\t')
    imm_df = pd.DataFrame(imm_reader)
    # Drop 1st column in order for the merge to function properly
    imm_df.drop(columns=['Sample_ID'], inplace = True)
    # Rename columns in order for the merge to function properly
    imm_df.rename(columns={"mut_peptide":"MT Epitope Seq",
                                    "Reference":"WT Epitope Seq",
                                    "peptide_variant_position":"Mutation Position",
                                    "TCGA_predict":"Immunogenicity_score"}, inplace=True)
    # Inner join dataFrames
    merged_df = pd.merge(pvacseq_df, imm_df)
    merged_df.to_csv(output, sep="\t", index=False)

if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(description='Calculate CSiN')

    parser.add_argument('--pvacseq_tsv', required=False, help='Input filtered MHC I pVACseq TSV file')
    parser.add_argument('--score_tsv', required=False, help='Input immunogenicity scores TSV files')
    parser.add_argument('--output', required=True, help='Path to output file')

    args = parser.parse_args()
    # Parse arguments
    pvacseq_tsv = args.pvacseq_tsv
    immun_tsv = args.score_tsv
    output = args.output

    append_score(pvacseq_tsv, immun_tsv, output)