#!/usr/bin/env python

"""
Parse the output of OptiType and HLA-HD and get HLA types (I and II) at 4 digit resolution

Requirements:
    * Python >= 2.7

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
    class_I = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-V"]
    class_II = []
    for hla in list(hlas):
        hla_gene = hla.split("*")[0]
        if hla_gene in class_I:
            continue
        else:
            class_II.append(hla)

    return class_II


def parse_hlahd(inFile, hlas=[]):
    hlas_tmp = []
    csv_reader = csv.reader(inFile, delimiter="\t")
    # inFile.readline()
    for line in csv_reader:
        if "HLA-" in line[1]:
            if len(line[1].split("-")[1].split(":")[0].split("*")[0]) > 1:
                hlas_tmp.append(line[1].split("-")[1].split(":")[0] + ":" + line[1].split(":")[1])
            else:
                hlas_tmp.append(line[1].split(":")[0] + ":" + line[1].split(":")[1])
        if "HLA-" in line[2]:
            if len(line[2].split("-")[1].split(":")[0].split("*")[0]) > 1:
                hlas_tmp.append(line[2].split("-")[1].split(":")[0] + ":" + line[2].split(":")[1])
            else:
                hlas_tmp.append(line[2].split(":")[0] + ":" + line[2].split(":")[1])

    for hla in filter_class_I(hlas_tmp):
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


def merge_hlas(hlas_DNA=[], hlas_RNA=[]):
    """
    From slack discussion on 20200425 (FF, GF, DR)
    We run Optitype RNA and WES in the main pipeline (only on tumor)
    We consider the WES class I HLA
    IFF the WES-based HLA are homo and are a subset of RNA-based HLA, then we consider also the second HLA alleles predicted from RNA
    These class I HLA are used for the somatic pipeline and also for the embedded NeoFuse
    """

    hlas = []
    hla_sets_RNA = {"HLA-A": [], "HLA-B": [], "HLA-C": []}
    hla_sets_DNA = {"HLA-A": [], "HLA-B": [], "HLA-C": []}

    for hla in hlas_RNA:
        hla_gene = hla.split("*")[0]
        hla_sets_RNA[hla_gene].append(hla)

    for hla in hlas_DNA:
        hla_gene = hla.split("*")[0]
        hla_sets_DNA[hla_gene].append(hla)

    for hla_gene in hla_sets_DNA.keys():
        if set(hla_sets_DNA[hla_gene]) == set(hla_sets_RNA[hla_gene]):
            hlas.append(hla_sets_DNA[hla_gene])
        elif hla_sets_DNA[hla_gene][0] == hla_sets_DNA[hla_gene][1]:
            if hla_sets_DNA[hla_gene][0] in hla_sets_RNA[hla_gene]:
                hlas.append(hla_sets_RNA[hla_gene])
        else:
            hlas.append(hla_sets_DNA[hla_gene])

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
    parser.add_argument("--opti_out", required=True, type=_file_read, help="output file produced by OptiType")
    parser.add_argument(
        "--opti_out_RNA",
        required=False,
        type=_file_read,
        default=False,
        help="output file produced by OptiType from RNAseq",
    )
    parser.add_argument("--hlahd_out", required=True, type=_file_read, help="output file produced by HLA-HD")
    parser.add_argument(
        "--custom", required=False, type=_file_read, default=False, help="User supplied 4 digit HLA type file"
    )
    parser.add_argument("--ref_hlas", required=True, type=_file_read, help="valid HLAs for pVACseq")

    args = parser.parse_args()

    optitype_dna_result = args.opti_out
    optitype_rna_result = args.opti_out_RNA
    hlahd_dna_result = args.hlahd_out
    reference_hlas = args.ref_hlas

    hla_array = []

    parse_opti(optitype_dna_result, hla_array)

    if args.opti_out_RNA is not False:
        hla_array_RNA = []
        parse_opti(optitype_rna_result, hla_array_RNA)
        if hla_array_RNA != []:
            hla_array = merge_hlas(hla_array, hla_array_RNA)

    parse_hlahd(hlahd_dna_result, hla_array)

    if args.custom is not False:
        read_custom_hla(args.custom, hla_array)

    valid_hlas = check_hlas(reference_hlas, hla_array)

    for hla in list(set(valid_hlas)):
        print(hla)
