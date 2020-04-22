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


def parse_hlahd(inFile, hlas=[]):
    csv_reader = csv.reader(inFile, delimiter="\t")
    inFile.readline()
    for line in csv_reader:
        if "HLA-" in line[1]:
            if len(line[1].split("-")[1].split(":")[0].split("*")[0]) > 1:
                hlas.append(line[1].split("-")[1].split(":")[0] + ":" + line[1].split(":")[1])
            else:
                hlas.append(line[1].split(":")[0] + ":" + line[1].split(":")[1])
        if "HLA-" in line[2]:
            if len(line[2].split("-")[1].split(":")[0].split("*")[0]) > 1:
                hlas.append(line[2].split("-")[1].split(":")[0] + ":" + line[2].split(":")[1])
            else:
                hlas.append(line[2].split(":")[0] + ":" + line[2].split(":")[1])
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
    parser.add_argument("--hlahd_out", required=True, type=_file_read, help="output file produced by HLA-HD")
    parser.add_argument(
        "--custom", required=False, type=_file_read, default=False, help="User supplied 4 digit HLA type file"
    )
    parser.add_argument("--ref_hlas", required=True, type=_file_read, help="valid HLAs for pVACseq")

    args = parser.parse_args()

    infile = args.opti_out
    infile2 = args.hlahd_out
    infile3 = args.ref_hlas

    hla_array = []

    parse_opti(infile, hla_array)
    parse_hlahd(infile2, hla_array)
    if args.custom is not False:
        read_custom_hla(args.custom, hla_array)

    valid_hlas = check_hlas(infile3, hla_array)

    for hla in list(set(valid_hlas)):
        print(hla)
