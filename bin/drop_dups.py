#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2
    * Pandas
    * NumPy

Copyright (c) 2020 Georgios Fotakis <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '3', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''

import argparse
import math
import os,sys
import pandas as pd
import numpy as np

if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(description='Calculate CSiN')

    parser.add_argument('--inFile', required=True, help='Input unfiltered pVACseq merged TSV files')
    parser.add_argument('--out', required=True, help='Input unfiltered pVACseq merged TSV files')

    args = parser.parse_args()

    inFile = args.inFile
    out = args.out

    inFile_reader = pd.read_csv(inFile, sep = '\t')
    inFile_df = pd.DataFrame(inFile_reader)
    no_dups_df = inFile_df.drop_duplicates()
    # no_dups_df = inFile_df.loc[(inFile_df['Chromosome'] == "chr2") & (inFile_df['Start'] == 43225325) & (inFile_df['Stop'] <= 43225325)]

    no_dups_df.to_csv(out, sep = "\t", index = False)