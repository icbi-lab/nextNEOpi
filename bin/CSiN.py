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
__version_info__ = (
    "0",
    "3",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""

import argparse
import math
import os, sys
import pandas as pd
import numpy as np


def convert_to_df(pvacseq_1_tsv, pvacseq_2_tsv):
    if not pvacseq_2_tsv:
        pvacseq_1_reader = pd.read_csv(pvacseq_1_tsv, sep="\t")
        merged_df = pd.DataFrame(pvacseq_1_reader)
    elif not pvacseq_1_tsv:
        pvacseq_2_reader = pd.read_csv(pvacseq_2_tsv, sep="\t")
        merged_df = pd.DataFrame(pvacseq_2_reader)
        # Rename columns in order for the filters to work
        merged_df.rename(columns={"NetMHCIIpan MT Percentile": "NetMHCpan MT Percentile"}, inplace=True)
    else:
        pvacseq_1_reader = pd.read_csv(pvacseq_1_tsv, sep="\t")
        pvacseq_1_df = pd.DataFrame(pvacseq_1_reader)
        pvacseq_2_reader = pd.read_csv(pvacseq_2_tsv, sep="\t")
        pvacseq_2_df = pd.DataFrame(pvacseq_2_reader)
        # Rename columns in order for the merge to work properly
        pvacseq_2_df.rename(
            columns={
                "NetMHCIIpan WT IC50 Score": "NetMHCpan WT IC50 Score",
                "NetMHCIIpan MT IC50 Score": "NetMHCpan MT IC50 Score",
                "NetMHCIIpan WT Percentile": "NetMHCpan WT Percentile",
                "NetMHCIIpan MT Percentile": "NetMHCpan MT Percentile",
            },
            inplace=True,
        )
        merged_df = pd.concat([pvacseq_1_df, pvacseq_2_df])
    return merged_df


def sub_csin(c, IC50_cutoff, xp_cutoff, filtered_df):
    if filtered_df.empty:
        return None
    else:
        sub_csins = []
        for rank in c:
            # Filter dataframe
            filtered_df_tmp = filtered_df
            filtered_df_tmp = filtered_df_tmp[filtered_df_tmp["NetMHCpan MT Percentile"] < rank]
            filtered_df_tmp = filtered_df_tmp[filtered_df_tmp["Best MT IC50 Score"] < IC50_cutoff]
            filtered_df_tmp = filtered_df_tmp[filtered_df_tmp["Gene Expression"] > xp_cutoff]
            # Get the VAF and mean VAF, then normilize
            vaf_mean = filtered_df_tmp["Tumor DNA VAF"].mean()
            filtered_df_tmp["Normalized_VAF"] = filtered_df_tmp["Tumor DNA VAF"].div(vaf_mean)
            # Get the load and mean load, then normilize
            temp = filtered_df_tmp.groupby(["Chromosome", "Start", "Stop"]).size().reset_index(name="vcf_load_count")
            vcf_load = temp["vcf_load_count"].mean()
            filtered_df_tmp["Normalized_load"] = filtered_df_tmp.groupby("Start")["Start"].transform("count").div(vcf_load)
            # Calculate "Sub-CSiN" and multiply by the a (=10) constant
            mult = filtered_df_tmp.Normalized_VAF * filtered_df_tmp.Normalized_load
            sub_csin = math.log(mult.mean()) * 10
            if math.isnan(sub_csin):
                sub_csin = 0.0
            else:
                pass
            sub_csins.append(sub_csin)
        return sub_csins


def csin(output, mhci_sub=None, mhcii_sub=None, merged_sub=None):
    # Claculate CSiN for each MHC class and for merged sub-CSiN DF
    if not mhci_sub:
        mhc_i, mhc_ii, merged = "NA", np.mean(mhcii_sub), "NA"
    elif not mhcii_sub:
        mhc_i, mhc_ii, merged = np.mean(mhci_sub), "NA", "NA"
    else:
        mhc_i, mhc_ii, merged = np.mean(mhci_sub), np.mean(mhcii_sub), np.mean(merged_sub)
    # Print final output
    with open(output, "w") as out:
        out.write("MHC I CSiN = %s\n\nMHC II CSiN = %s\n\nTotal CSiN = %s" % (mhc_i, mhc_ii, merged))


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Calculate CSiN")

    parser.add_argument("--MHCI_tsv", required=False, help="Input unfiltered pVACseq merged TSV files")
    parser.add_argument("--MHCII_tsv", required=False, help="Input unfiltered pVACseq merged TSV files")
    parser.add_argument(
        "--rank",
        required=True,
        type=float,
        nargs="+",
        help="The maximum cutoff rank value for an epitope to be considered as an HLA binder",
    )
    parser.add_argument(
        "--ic50",
        required=True,
        type=float,
        help="The maximum cutoff value for an epitope to be considered as an HLA binder",
    )
    parser.add_argument(
        "--gene_exp",
        required=True,
        type=float,
        help="The minimum expression cutoff value for a gene to be considered expressed",
    )
    parser.add_argument("--output", required=True, help="Path to the output directory")

    args = parser.parse_args()
    # Parse arguments
    MHCI_tsv = args.MHCI_tsv
    MHCII_tsv = args.MHCII_tsv
    ranks = args.rank
    ic50_cutoff = args.ic50
    xp_cutoff = args.gene_exp
    output = args.output

    # Create the DataFrames (depending on user input)
    if not MHCII_tsv:
        mhci_df = convert_to_df(MHCI_tsv, MHCII_tsv)
        mhcii_df = pd.DataFrame()
        merged_df = pd.DataFrame()
    elif not MHCI_tsv:
        mhci_df = pd.DataFrame()
        mhcii_df = convert_to_df(MHCI_tsv, MHCII_tsv)
        merged_df = pd.DataFrame()
    else:
        mhci = None
        mhcii = None
        mhci_df = convert_to_df(MHCI_tsv, mhcii)
        mhcii_df = convert_to_df(mhci, MHCII_tsv)
        merged_df = convert_to_df(MHCI_tsv, MHCII_tsv)

    # Calculate sub-CSiN
    filtered_mhci = sub_csin(ranks, ic50_cutoff, xp_cutoff, mhci_df)
    filtered_mhcii = sub_csin(ranks, ic50_cutoff, xp_cutoff, mhcii_df)
    filtered_merged = sub_csin(ranks, ic50_cutoff, xp_cutoff, merged_df)

    # Calculate CSiN
    csin(output, filtered_mhci, filtered_mhcii, filtered_merged)
