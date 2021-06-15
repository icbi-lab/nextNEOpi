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
import os, sys
import pandas as pd
from functools import reduce


def create_neofuse_consensus(in_files, out_file1, out_file2):
    df_arr = []
    new_df_arr = []
    for in_file in in_files:
        pd_reader = pd.read_csv(in_file, sep="\t")
        df = pd.DataFrame(pd_reader)
        df_arr.append(df)
    NeoFuse_keys = [
        "Fusion",
        "Gene1",
        "Gene2",
        "HLA_Type",
        "Fusion_Peptide",
        # "IC50",
        "Rank",
        "Event_Type",
        "Stop_Codon",
        "Confidence",
    ]
    Neofuse_drop = [
        "Gene1_TPM_x",
        "Gene2_TPM_x",
        "Avg_TPM_x",
        "HLA_TPM_x",
        "Gene1_TPM_y",
        "Gene2_TPM_y",
        "Avg_TPM_y",
        "HLA_TPM_y",
        "Gene1_TPM_x_avg",
        "Gene2_TPM_x_avg",
        "Gene1_TPM_y_avg",
        "Gene2_TPM_y_avg",
        "Avg_TPM_xx",
        "Avg_TPM_yy",
        "HLA_TPM_xx",
        "HLA_TPM_yy",
    ]

    intersection = reduce(lambda left, right: pd.merge(left, right, on=NeoFuse_keys), df_arr).drop_duplicates()
    # Calculate avg TPM for Gene1
    intersection["Gene1_TPM_x_avg"] = intersection["Gene1_TPM_x"].mean(axis=1)
    intersection["Gene1_TPM_y_avg"] = intersection["Gene1_TPM_y"].mean(axis=1)
    intersection["Gene1_avg_TPM"] = (intersection["Gene1_TPM_x_avg"] + intersection["Gene1_TPM_y_avg"]).div(2).round(2)
    # Calculate avg TPM for Gene2
    intersection["Gene2_TPM_x_avg"] = intersection["Gene2_TPM_x"].mean(axis=1)
    intersection["Gene2_TPM_y_avg"] = intersection["Gene2_TPM_y"].mean(axis=1)
    intersection["Gene2_avg_TPM"] = (intersection["Gene2_TPM_x_avg"] + intersection["Gene2_TPM_y_avg"]).div(2).round(2)
    # Calculate avg TPM for the fusion gene
    intersection["Avg_TPM_xx"] = intersection["Avg_TPM_x"].mean(axis=1)
    intersection["Avg_TPM_yy"] = intersection["Avg_TPM_y"].mean(axis=1)
    intersection["Avg_TPM"] = (intersection["Avg_TPM_xx"] + intersection["Avg_TPM_yy"]).div(2).round(2)
    # Calculate avg TPM for the HLA types
    intersection["HLA_TPM_xx"] = intersection["HLA_TPM_x"].mean(axis=1)
    intersection["HLA_TPM_yy"] = intersection["HLA_TPM_y"].mean(axis=1)
    intersection["HLA_avg_TPM"] = (intersection["Avg_TPM_xx"] + intersection["Avg_TPM_yy"]).div(2).round(2)
    # Drop unwanted columns and print to final output file
    intersection = intersection.drop(Neofuse_drop, axis=1)
    intersection.to_csv(out_file1, index=False, sep="\t")

    union = pd.concat(df_arr, ignore_index=True).drop_duplicates().reset_index(drop=True)
    union.to_csv(out_file2, index=False, sep="\t")

    return 0


def create_pvacseq_consensus(in_files, out_file1, out_file2, in_type):

    df_arr = []
    new_df_arr = []
    files_with_data = []
    for in_file in in_files:
        pd_reader = pd.read_csv(in_file, sep="\t")
        df = pd.DataFrame(pd_reader)
        if len(df) > 0:
            df_arr.append(df)
            files_with_data.append(in_file)

    if len(df_arr) > 1:
        print("Making intersect and union of: " + " ".join(files_with_data))
        pvacseq_keys = [
            "Chromosome",
            "Start",
            "Stop",
            "Reference",
            "Variant",
            "Transcript",
            "Transcript Support Level",
            "Ensembl Gene ID",
            "Variant Type",
            "Mutation",
            "Protein Position",
            "Gene Name",
            # "HGVSc",
            # "HGVSp",
            "HLA Allele",
            "Peptide Length",
            "Sub-peptide Position",
            "Mutation Position",
            "MT Epitope Seq",
            "WT Epitope Seq",
            "Best MT Score Method",
            "Best MT Score",
            # "Corresponding WT Score",
            # "Corresponding Fold Change",
            "Best MT Percentile Method",
            "Best MT Percentile",
            # "Corresponding WT Percentile",
            # "Transcript Expression",
            "Median MT Score",
            # "Median WT Score",
            # "Median Fold Change",
            "Median MT Percentile",
            "Median WT Percentile",
            # "cterm_7mer_gravy_score",
            # "max_7mer_gravy_score",
            # "difficult_n_terminal_residue",
            # "c_terminal_cysteine",
            # "c_terminal_proline",
            # "cysteine_count",
            # "n_terminal_asparagine",
            # "asparagine_proline_bond_count"
        ]
        pvacseq_drop = [
            "Tumor DNA Depth_x",
            "Tumor DNA Depth_y",
            "Tumor DNA VAF_x",
            "Tumor DNA VAF_y",
            "Tumor RNA Depth_x",
            "Tumor RNA Depth_y",
            "Tumor RNA VAF_x",
            "Tumor RNA VAF_y",
            "Normal Depth_x",
            "Normal Depth_y",
            "Normal VAF_x",
            "Normal VAF_y",
            "Gene Expression_x",
            "Gene Expression_y",
            "Index_x",
            "Index_y",
            "Tumor DNA Depth_x_avg",
            "Tumor DNA Depth_y_avg",
            "Tumor DNA VAF_x_avg",
            "Tumor DNA VAF_y_avg",
            "Tumor RNA Depth_x_avg",
            "Tumor RNA Depth_y_avg",
            "Tumor RNA VAF_x_avg",
            "Tumor RNA VAF_y_avg",
            "Normal Depth_x_avg",
            "Normal Depth_y_avg",
            "Normal VAF_x_avg",
            "Normal VAF_y_avg",
            "Gene Expression_x_avg",
            "Gene Expression_y_avg",
        ]
        if in_type == "mhc_i":
            mhc_i_cols = [
                "NetMHCpan WT Score",
                "NetMHCpan MT Score",
                "NetMHCpan WT Percentile",
                "NetMHCpan MT Percentile",
            ]
            pvacseq_keys = pvacseq_keys + mhc_i_cols
        elif in_type == "mhc_ii":
            mhc_ii_cols = [
                "NetMHCIIpan WT Score",
                "NetMHCIIpan MT Score",
                "NetMHCIIpan WT Percentile",
                "NetMHCIIpan MT Percentile",
            ]
            pvacseq_keys = pvacseq_keys + mhc_ii_cols
        intersection = reduce(lambda left, right: pd.merge(left, right, on=pvacseq_keys), df_arr).drop_duplicates()
        # Calculate avg Tumor DNA Depth
        intersection["Tumor DNA Depth_x_avg"] = intersection["Tumor DNA Depth_x"].mean(axis=1)
        intersection["Tumor DNA Depth_y_avg"] = intersection["Tumor DNA Depth_y"].mean(axis=1)
        intersection["Tumor DNA Depth_avg"] = (
            (intersection["Tumor DNA Depth_x_avg"] + intersection["Tumor DNA Depth_y_avg"]).div(2).round(2)
        )
        # Calculate avg Tumor DNA VAF
        intersection["Tumor DNA VAF_x_avg"] = intersection["Tumor DNA VAF_x"].mean(axis=1)
        intersection["Tumor DNA VAF_y_avg"] = intersection["Tumor DNA VAF_y"].mean(axis=1)
        intersection["Tumor DNA VAF_avg"] = (
            (intersection["Tumor DNA VAF_x_avg"] + intersection["Tumor DNA VAF_y_avg"]).div(2).round(2)
        )
        # Calculate avg Tumor RNA Depth
        intersection["Tumor RNA Depth_x_avg"] = intersection["Tumor RNA Depth_x"].mean(axis=1)
        intersection["Tumor RNA Depth_y_avg"] = intersection["Tumor RNA Depth_y"].mean(axis=1)
        intersection["Tumor RNA Depth_avg"] = (
            (intersection["Tumor RNA Depth_x_avg"] + intersection["Tumor RNA Depth_y_avg"]).div(2).round(2)
        )
        # Calculate avg Tumor RNA VAF
        intersection["Tumor RNA VAF_x_avg"] = intersection["Tumor RNA VAF_x"].mean(axis=1)
        intersection["Tumor RNA VAF_y_avg"] = intersection["Tumor RNA VAF_y"].mean(axis=1)
        intersection["Tumor RNA VAF_avg"] = (
            (intersection["Tumor RNA VAF_x_avg"] + intersection["Tumor RNA VAF_y_avg"]).div(2).round(2)
        )
        # Calculate avg Normal Depth
        intersection["Normal Depth_x_avg"] = intersection["Normal Depth_x"].mean(axis=1)
        intersection["Normal Depth_y_avg"] = intersection["Normal Depth_y"].mean(axis=1)
        intersection["Normal Depth_avg"] = (
            (intersection["Normal Depth_x_avg"] + intersection["Normal Depth_y_avg"]).div(2).round(2)
        )
        # Calculate avg Normal VAF
        intersection["Normal VAF_x_avg"] = intersection["Normal VAF_x"].mean(axis=1)
        intersection["Normal VAF_y_avg"] = intersection["Normal VAF_y"].mean(axis=1)
        intersection["Normal VAF_avg"] = (
            (intersection["Normal VAF_x_avg"] + intersection["Normal VAF_y_avg"]).div(2).round(2)
        )
        # Calculate avg Gene Expression
        intersection["Gene Expression_x_avg"] = intersection["Gene Expression_x"].mean(axis=1)
        intersection["Gene Expression_y_avg"] = intersection["Gene Expression_y"].mean(axis=1)
        intersection["Gene Expression_avg"] = (
            (intersection["Gene Expression_x_avg"] + intersection["Gene Expression_y_avg"]).div(2).round(2)
        )
        # Drop unwanted columns and print to final intersection output file
        intersection = intersection.drop(pvacseq_drop, axis=1)
        intersection.to_csv(out_file1, index=False, sep="\t")
        # Drop unwanted columns and print to final uninon output file
        union = pd.concat(df_arr, ignore_index=True).drop_duplicates().reset_index(drop=True)
        union.to_csv(out_file2, index=False, sep="\t")
    else:
        if len(files_with_data) == 1:
            print("Less than two files with data content: can not make intersect or union")
            print("File with data rows: " + files_with_data[0])
        else:
            print("No data rows in inputfiles")


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Create consesus list from final pVACseq output")

    parser.add_argument("--input_files", required=True, help="Input pVACseq merged TSV files")

    parser.add_argument("--output", required=True, help="Output file")

    parser.add_argument("--pvacseq", nargs="?", const=True, default=False, help="Activate column parsing mode.")

    args = parser.parse_args()

    # Parse arguments
    files = args.input_files.split(",")
    out_file_name = args.output
    pvacseq = args.pvacseq
    # Function calls
    if not pvacseq:
        output1 = out_file_name + "_NeoFuse_intersection.tsv"
        output2 = out_file_name + "_NeoFuse_union.tsv"
        create_neofuse_consensus(files, output1, output2)
    elif pvacseq == "mhc_i":
        output1 = out_file_name + "_MHCI_intersection.tsv"
        output2 = out_file_name + "_MHCI_union.tsv"
        create_pvacseq_consensus(files, output1, output2, pvacseq)
    elif pvacseq == "mhc_ii":
        output1 = out_file_name + "_MHCII_intersection.tsv"
        output2 = out_file_name + "_MHCII_union.tsv"
        create_pvacseq_consensus(files, output1, output2, pvacseq)
