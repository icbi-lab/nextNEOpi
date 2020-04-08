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

def merge(pvacseq_1_tsv, pvacseq_2_tsv = None):
    if not pvacseq_2_tsv:
        pvacseq_1_reader = pd.read_csv(pvacseq_1_tsv, sep = '\t')
        merged_df = pd.DataFrame(pvacseq_1_reader)
    elif not pvacseq_1_tsv:
        pvacseq_2_reader = pd.read_csv(pvacseq_2_tsv, sep = '\t')
        merged_df = pd.DataFrame(pvacseq_2_reader)
        merged_df.rename(columns={"NetMHCIIpan WT Score":"NetMHCpan WT Score",
                                    "NetMHCIIpan MT Score":"NetMHCpan MT Score",
                                    "NetMHCIIpan WT Percentile":"NetMHCpan WT Percentile",
                                    "NetMHCIIpan MT Percentile":"NetMHCpan MT Percentile"}, inplace=True)
    else:
        pvacseq_1_reader = pd.read_csv(pvacseq_1_tsv, sep = '\t')
        pvacseq_1_df = pd.DataFrame(pvacseq_1_reader)
        pvacseq_2_reader = pd.read_csv(pvacseq_2_tsv, sep = '\t')
        pvacseq_2_df = pd.DataFrame(pvacseq_2_reader)
        pvacseq_2_df.rename(columns={"NetMHCIIpan WT Score":"NetMHCpan WT Score",
                                    "NetMHCIIpan MT Score":"NetMHCpan MT Score",
                                    "NetMHCIIpan WT Percentile":"NetMHCpan WT Percentile",
                                    "NetMHCIIpan MT Percentile":"NetMHCpan MT Percentile"}, inplace=True)
        merged_df = pd.concat([pvacseq_1_df, pvacseq_2_df])
    # print(merged_tsv)
    return merged_df



def filter_df(pvacseq_df, c, IC50_cutoff, xp_cutoff, output):
    # import pVACseq tsv as pandas DataFrame
    # pvacseq_reader = pd.read_csv(pvacseq_tsv, sep = '\t')
    # pvacseq_df = pd.DataFrame(pvacseq_reader)
    pvacseq_df = pvacseq_df[pvacseq_df['NetMHCpan MT Percentile'] < c]
    pvacseq_df = pvacseq_df[pvacseq_df['Best MT Score'] < IC50_cutoff]
    pvacseq_df = pvacseq_df[pvacseq_df['Gene Expression'] > xp_cutoff]

    # DEBUG: print filtered data frames - for dev purposes only - TODO: delete on stable version
    print("Resulting DataFrame after filtering for [c < %s, gene expression > %s and IC50 < %s]:\n\n%s\n\nWriting to file: %s" % 
            (c, xp_cutoff, IC50_cutoff, pvacseq_df, os.path.realpath(output)))
    # DEBUG: create temp output files to check the results - for dev purposes only - TODO: delete on stable version
    with open(output, "w") as out:
        pvacseq_df.to_csv(out, sep = "\t", index = False)
    
    return pvacseq_df

def csin(filtered_df):
    # Get the VAF and mean VAF, then normilize
    vaf_mean = filtered_df['Tumor DNA VAF'].mean()
    filtered_df['Normalized_VAF'] = filtered_df['Tumor DNA VAF'].div(vaf_mean)
    # Get the load and mean load, then normilize
    temp = filtered_df.groupby(['Start']).size().reset_index(name='vcf_load_count')
    vcf_load = temp['vcf_load_count'].mean()
    filtered_df['Normalized_load'] = filtered_df.groupby('Start')['Start'].transform('count').div(vcf_load)
    # Calculate "Sub-CSiN" and multiply by the a (=10) constant
    mult = filtered_df.Normalized_VAF * filtered_df.Normalized_load
    sub_csin = math.log(mult.mean())*10
    # DEBUG: print CSiN value - for dev purposes only - TODO: delete on stable version
    print("Resulting Sub-CSiN: %s\n" % sub_csin)

    return sub_csin

if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(description='Calculate CSiN')

    parser.add_argument('--MHCI_tsv', required=False, help='Input unfiltered pVACseq merged TSV files')
    parser.add_argument('--MHCII_tsv', required=False, help='Input unfiltered pVACseq merged TSV files')
    parser.add_argument('--rank', required=True, type = float, nargs='+', help='The maximum cutoff rank value for an epitope to be considered as an HLA binder')
    parser.add_argument('--ic50', required=True, type = float, help='The maximum cutoff value for an epitope to be considered as an HLA binder')
    parser.add_argument('--gene_exp', required=True, type = float, help='The minimum expression cutoff value for a gene to be considered expressed')
    parser.add_argument('--output', required=True, help='Path to the output directory')

    args = parser.parse_args()

    MHCI_tsv = args.MHCI_tsv
    MHCII_tsv = args.MHCII_tsv
    ranks = args.rank
    ic50_cutoff = args.ic50
    xp_cutoff = args.gene_exp
    output = args.output


    # Create merged DataFrame
    merged_df = merge(MHCI_tsv, MHCII_tsv)

    # Store sub-CSiN in array
    sub_csin = []
    for c in ranks:
        filtered_df = filter_df(merged_df, c, ic50_cutoff, xp_cutoff, output+"CSiN_%s.tsv" % c)
        sub_csin.append(csin(filtered_df))

    # Calculate and print CSiN
    print("CSiN = %s" % np.mean(sub_csin))

