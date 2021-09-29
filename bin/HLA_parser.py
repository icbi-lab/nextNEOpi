#!/usr/bin/env python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2021 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

import argparse
import os
import sys
import csv


def flatten(l):
    results = []
    for rec in l:
        if type(rec) == list:
            results += rec
            results = flatten(results)
        else:
            results.append(rec)
    return results


def filter_class_I(hlas=[]):
    ref_I = [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "HLA-F",
        "HLA-G",
        "HLA-H",
        "HLA-J",
        "HLA-K",
        "HLA-L",
        "HLA-V",
    ]
    class_II = []
    for hla in list(hlas):
        hla_gene = hla.split("*")[0]
        if hla_gene in ref_I:
            continue
        else:
            class_II.append(hla)
    return class_II

def filter_class_II(hlas=[]):
    ref_I = [
        "HLA-A",
        "HLA-B",
        "HLA-C"
    ]
    class_I = []
    for hla in list(hlas):
        hla_gene = hla.split("*")[0]
        if hla_gene in ref_I:
            class_I.append(hla)
        else:
            continue
    return class_I

def parse_hlahd(inFile,
                hla_class,
                hlas=[]):
    hlas_tmp = []
    csv_reader = csv.reader(inFile, delimiter="\t")
    # inFile.readline()
    for line in csv_reader:
        if "HLA-" in line[1]:
            if len(line[1].split("-")[1].split(":")[0].split("*")[0]) > 1:
                hlas_tmp.append(
                    line[1].split("-")[1].split(":")[0] + ":" + line[1].split(":")[1]
                    # line[1].split(":")[0] + ":" + line[1].split(":")[1]
                )
            else:
                hlas_tmp.append(line[1].split(":")[0] + ":" + line[1].split(":")[1])
        if "HLA-" in line[2]:
            if len(line[2].split("-")[1].split(":")[0].split("*")[0]) > 1:
                hlas_tmp.append(
                    line[2].split("-")[1].split(":")[0] + ":" + line[2].split(":")[1]
                    # line[2].split(":")[0] + ":" + line[2].split(":")[1]
                )
            else:
                hlas_tmp.append(line[2].split(":")[0] + ":" + line[2].split(":")[1])
    inFile.seek(0)

    if hla_class == "classII":
        for hla in filter_class_I(hlas_tmp):
            hlas.append(hla)
    elif hla_class == "classI":
        for hla in filter_class_II(hlas_tmp):
            hlas.append(hla)

    return hlas


def parse_opti(inFile, hlas=[]):
    csv_reader = csv.reader(inFile, delimiter="\t")
    inFile.readline()
    for line in csv_reader:
        if "*" in line[1]:
            hlas.append("HLA-" + line[1])
        if "*" in line[2]:
            hlas.append("HLA-" + line[2])
        if "*" in line[3]:
            hlas.append("HLA-" + line[3])
        if "*" in line[4]:
            hlas.append("HLA-" + line[4])
        if "*" in line[5]:
            hlas.append("HLA-" + line[5])
        if "*" in line[6]:
            hlas.append("HLA-" + line[6])
    return hlas


def merge_hlasI(hlas_DNA=[], hlas_RNA=[]):
    """
    From slack discussion on 20200425 (FF, GF, DR).
    We run Optitype RNA and WES in the main pipeline (only on tumor).
    We consider the WES class I HLA.
    IF the WES-based HLA are homo and are a subset of RNA-based HLA, then we consider also the second HLA alleles predicted from RNA.
    These class I HLA are used for the somatic pipeline and also for the embedded NeoFuse.
    """

    hlas = []

    hlaI_sets_RNA = {"HLA-A": [], "HLA-B": [], "HLA-C": []}

    hlaI_sets_DNA = {"HLA-A": [], "HLA-B": [], "HLA-C": []}

    for hla in hlas_RNA:
        hla_gene = hla.split("*")[0]
        hlaI_sets_RNA[hla_gene].append(hla)

    for hla in hlas_DNA:
        hla_gene = hla.split("*")[0]
        hlaI_sets_DNA[hla_gene].append(hla)

    for hla_gene in hlaI_sets_DNA.keys():
        if set(hlaI_sets_DNA[hla_gene]) == set(hlaI_sets_RNA[hla_gene]):
            hlas.append(hlaI_sets_DNA[hla_gene])
        elif len(hlaI_sets_DNA[hla_gene]) == 0:
            if len(hlaI_sets_RNA[hla_gene]) != 0:
                hlas.append(hlaI_sets_RNA[hla_gene])
        elif (
            len(hlaI_sets_DNA[hla_gene]) == 1
            or hlaI_sets_DNA[hla_gene][0] == hlaI_sets_DNA[hla_gene][1]
        ):
            if hlaI_sets_DNA[hla_gene][0] in hlaI_sets_RNA[hla_gene]:
                hlas.append(hlaI_sets_RNA[hla_gene])
            else:
                hlas.append(hlaI_sets_DNA[hla_gene])
            #     hlas.append(hlaI_sets_RNA[hla_gene])
        else:
            hlas.append(hlaI_sets_DNA[hla_gene])

    return flatten(hlas)


def merge_hlasII(hlas_DNA=[], hlas_RNA=[]):
    """
    From slack discussion on 20200425 (FF, GF, DR).
    We run HLA-HD RNA and WES in the main pipeline (only on tumor).
    We consider the WES class II HLA.
    IF the WES-based HLA II are homo and are a subset of RNA-based HLA, then we consider also the second HLA II alleles predicted from RNA.
    These class II HLA are used for the somatic pipeline and also for the embedded NeoFuse.
    """

    hlasII = []

    hlaII_sets_RNA = {
        "DRB1": [],
        "DQA1": [],
        "DQB1": [],
        "DPA1": [],
        "DPB1": [],
        "DMA": [],
        "DMB": [],
        "DOA": [],
        "DOB": [],
        "DRA": [],
        "DRB2": [],
        "DRB3": [],
        "DRB4": [],
        "DRB5": [],
        "DRB6": [],
        "DRB7": [],
        "DRB8": [],
        "DRB9": [],
    }

    hlaII_sets_DNA = {
        "DRB1": [],
        "DQA1": [],
        "DQB1": [],
        "DPA1": [],
        "DPB1": [],
        "DMA": [],
        "DMB": [],
        "DOA": [],
        "DOB": [],
        "DRA": [],
        "DRB2": [],
        "DRB3": [],
        "DRB4": [],
        "DRB5": [],
        "DRB6": [],
        "DRB7": [],
        "DRB8": [],
        "DRB9": [],
    }

    for hla in hlas_RNA:
        hla_gene = hla.split("*")[0]
        hlaII_sets_RNA[hla_gene].append(hla)

    for hla in hlas_DNA:
        hla_gene = hla.split("*")[0]
        hlaII_sets_DNA[hla_gene].append(hla)

    for hla_gene in hlaII_sets_DNA.keys():
        if set(hlaII_sets_DNA[hla_gene]) == set(hlaII_sets_RNA[hla_gene]):
            hlasII.append(hlaII_sets_DNA[hla_gene])
        elif len(hlaII_sets_DNA[hla_gene]) == 0:
            if len(hlaII_sets_RNA[hla_gene]) != 0:
                hlasII.append(hlaII_sets_RNA[hla_gene])
        elif (
            len(hlaII_sets_DNA[hla_gene]) == 1
            or hlaII_sets_DNA[hla_gene][0] == hlaII_sets_DNA[hla_gene][1]
        ):
            if hlaII_sets_DNA[hla_gene][0] in hlaII_sets_RNA[hla_gene]:
                hlasII.append(hlaII_sets_RNA[hla_gene])
            else:
                hlasII.append(hlaII_sets_DNA[hla_gene])
            #     hlasII.append(hlaII_sets_RNA[hla_gene])
        else:
            hlasII.append(hlaII_sets_DNA[hla_gene])

    return flatten(hlasII)


def forced_calls(hlas_I=[], hlas_II=[]):
    hlas = []

    hlaI_sets = {"HLA-A": [], "HLA-B": [], "HLA-C": []}

    hlaII_sets = {
        "DRB1": [],
        "DQA1": [],
        "DQB1": [],
        "DPA1": [],
        "DPB1": [],
        "DMA": [],
        "DMB": [],
        "DOA": [],
        "DOB": [],
        "DRA": [],
        "DRB2": [],
        "DRB3": [],
        "DRB4": [],
        "DRB5": [],
        "DRB6": [],
        "DRB7": [],
        "DRB8": [],
        "DRB9": [],
    }

    for hlaI in hlas_I:
        hla_gene = hlaI.split("*")[0]
        hlaI_sets[hla_gene].append(hlaI)

    for hlaII in hlas_II:
        hla_gene = hlaII.split("*")[0]
        hlaII_sets[hla_gene].append(hlaII)

    for hlaI_gene in hlaI_sets.keys():
        hlas.append(hlaI_sets[hlaI_gene])

    for hlaII_gene in hlaII_sets.keys():
        hlas.append(hlaII_sets[hlaII_gene])

    return flatten(hlas)


def read_custom_hla(inFile, hlas=[]):
    for line in inFile:
        line = line.strip()
        hlas.append(line)
    return hlas


def check_hlas(inFile, hlas=[]):
    known_hlas = []
    for line in inFile:
        line = line.strip()
        known_hlas.append(line)

    valid_hlas = []
    for hla in hlas:
        if hla in known_hlas:
            valid_hlas.append(hla)

    return valid_hlas


if __name__ == "__main__":

    def _file_read(fname, mode="rt"):
        """Returns an open file handle if the given filename exists."""
        if not os.path.exists(fname):
            parser.error("File '{0}' not found.".format(fname))
        _, file_extension = os.path.splitext(fname)
        if file_extension == ".gz":
            mode = "rb"
        return open(fname, mode)

    parser = argparse.ArgumentParser(
        description="Parse the output files of OptiType and HLA-HD for all HLA types and alleles"
    )
    parser.add_argument(
        "--opti_out",
        required=False,
        type=_file_read,
        help="output file produced by OptiType from WES",
    )
    parser.add_argument(
        "--opti_out_RNA",
        required=False,
        type=_file_read,
        default=False,
        help="output file produced by OptiType from RNAseq",
    )
    parser.add_argument(
        "--hlahd_out",
        required=False,
        type=_file_read,
        help="output file produced by HLA-HD from WES",
    )
    parser.add_argument(
        "--hlahd_out_RNA",
        required=False,
        type=_file_read,
        help="output file produced by HLA-HD from RNAseq",
    )
    parser.add_argument(
        "--custom",
        required=False,
        type=_file_read,
        default=False,
        help="User supplied 4 digit HLA type file",
    )
    parser.add_argument(
        "--ref_hlas", required=True, type=_file_read, help="valid HLAs for pVACseq"
    )
    parser.add_argument(
        "--force_DNA",
        required=False,
        action="store_true",
        help="return HLA types only from WES data",
    )
    parser.add_argument(
        "--force_RNA",
        required=False,
        action="store_true",
        help="return HLA types only from WES data",
    )

    args = parser.parse_args()

    # Get input
    optitype_dna_result = args.opti_out
    optitype_rna_result = args.opti_out_RNA
    hlahd_dna_result = args.hlahd_out
    hlahd_rna_result = args.hlahd_out_RNA
    reference_hlas = args.ref_hlas

    hlaI_array = []
    hlaII_array = []
    hla_array = []

    # Check if user defined only custom list
    if not optitype_dna_result and not hlahd_dna_result and args.custom:
        read_custom_hla(args.custom, hla_array)
    else:
        if args.force_DNA:
            if not optitype_dna_result:
                parse_hlahd(hlahd_dna_result, "classI", hlaI_array)
                parse_hlahd(hlahd_dna_result, "classII", hlaII_array)
                hla_array = forced_calls(hlaI_array, hlaII_array)
            else:
                parse_opti(optitype_dna_result, hlaI_array)
                parse_hlahd(hlahd_dna_result, "classII", hlaII_array)
                hla_array = forced_calls(hlaI_array, hlaII_array)
        elif args.force_RNA:
            if not optitype_rna_result:
                parse_hlahd(hlahd_rna_result, "classI" ,hlaI_array)
                parse_hlahd(hlahd_rna_result, "classII", hlaII_array)
                hla_array = forced_calls(hlaI_array, hlaII_array)
            else:
                parse_opti(optitype_rna_result, hlaI_array)
                parse_hlahd(hlahd_rna_result, "classII", hlaII_array)
                hla_array = forced_calls(hlaI_array, hlaII_array)
        else:
            if args.opti_out:
                parse_opti(optitype_dna_result, hlaI_array)
            else:
                parse_hlahd(hlahd_dna_result, "classI", hlaI_array)

            # Check if OptiType RNAseq data are available
            if args.opti_out_RNA:
                hla_array_opti_RNA = []

                parse_opti(optitype_rna_result, hla_array_opti_RNA)

                if hla_array_opti_RNA != []:
                    hlaI_array = merge_hlasI(hlaI_array, hla_array_opti_RNA)

            if args.hlahd_out:
                parse_hlahd(hlahd_dna_result, "classII", hlaII_array)

            # Check if HLA-HD RNAseq data are available
            if args.hlahd_out_RNA:
                hla_array_hlahd_RNA = []

                parse_hlahd(hlahd_rna_result, hla_array_hlahd_RNA)

                if hla_array_hlahd_RNA != []:
                    hlaII_array = merge_hlasII(hlaII_array, hla_array_hlahd_RNA)

            # Include custom list to final results if available
            if args.custom:
                read_custom_hla(args.custom, hla_array)

            # Merge HLA I and II arrays
            hla_array = hlaI_array + hlaII_array

    # Check for valid HLA types and return final result
    valid_hlas = check_hlas(reference_hlas, hla_array)

    for hla in list(set(valid_hlas)):
        print(hla)

