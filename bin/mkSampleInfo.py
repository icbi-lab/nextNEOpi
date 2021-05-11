#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.7
    * Pysam

Copyright (c) 2021 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>


"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""

import os
import sys
import argparse


def parse_csin(csin_fh, csin_info):
    for line in csin_fh:
        if line.find("MHC I CSiN") != -1:
            _, csin_v = line.split(" = ")
            csin_info["MHCI"] = round(float(csin_v.strip()), 3)
        if line.find("MHC II CSiN") != -1:
            _, csin_v = line.split(" = ")
            csin_info["MHCII"] = round(float(csin_v.strip()), 3)
        if line.find("Total CSiN") != -1:
            _, csin_v = line.split(" = ")
            csin_info["combined"] = round(float(csin_v.strip()), 3)

    return csin_info


def parse_tmb(tmb_fh, tmb_info, tmb_type):
    for line in tmb_fh:
        if line.find("Coverage") != -1:
            _, v = line.split("\t")
            if tmb_type == "all":
                tmb_info["cov_genome"] = v.strip()
            if tmb_type == "coding":
                tmb_info["cov_coding"] = v.strip()
        if line.find("Variants") != -1:
            _, v = line.split("\t")
            if tmb_type == "all":
                tmb_info["variants_tot"] = v.strip()
            if tmb_type == "coding":
                tmb_info["variants_coding"] = v.strip()
        if line.find("Mutational load (") != -1:
            _, v = line.split("\t")
            if tmb_type == "all":
                tmb_info["TMB"] = round(float(v.strip()), 3)
            if tmb_type == "coding":
                tmb_info["TMB_coding"] = round(float(v.strip()), 3)
        if line.find("Mutational load clonal") != -1:
            _, v = line.split("\t")
            if tmb_type == "all":
                tmb_info["TMB_clonal"] = round(float(v.strip()), 3)
            if tmb_type == "coding":
                tmb_info["TMB_clonal_coding"] = round(float(v.strip()), 3)


    return tmb_info


def write_output(out_fh, tmb_info, csin_info, sample_name):
    header_fields = [
        "SampleID",
        "TMB",
        "TMB_clonal",
        "TMB_coding",
        "TMB_clonal_coding",
        "variants_total",
        "variants_coding",
        "coverage_genome",
        "coverage_coding",
        "CSiN_MHC_I",
        "CSiN_MHC_II",
        "CSiN_combined",
    ]
    data_fields = [
        sample_name,
        tmb_info["TMB"],
        tmb_info["TMB_clonal"],
        tmb_info["TMB_coding"],
        tmb_info["TMB_clonal_coding"],
        tmb_info["variants_tot"],
        tmb_info["variants_coding"],
        tmb_info["cov_genome"],
        tmb_info["cov_coding"],
        csin_info["MHCI"],
        csin_info["MHCII"],
        csin_info["combined"],
    ]

    out_fh.write("\t".join(header_fields) + "\n")
    out_fh.write("\t".join(map(str, data_fields)) + "\n")


def _file_write(fname):
    """Returns an open file handle if the given filename exists."""
    return open(fname, "w")


def _file_read(fname):
    """Returns an open file handle if the given filename exists."""
    return open(fname, "r")


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Compile sample info sheet")
    parser.add_argument(
        "--tmb",
        required=True,
        type=_file_read,
        help="TMB file",
    )
    parser.add_argument(
        "--tmb_coding",
        required=True,
        type=_file_read,
        help="TMB coding file",
    )
    parser.add_argument(
        "--csin",
        required=True,
        type=_file_read,
        help="CSiN file",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=_file_write,
        help="Output file",
    )
    parser.add_argument(
        "--sample_name",
        required=True,
        type=str,
        help="Sample name",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args()

    tmb = args.tmb
    tmb_coding = args.tmb_coding
    csin = args.csin
    out = args.out
    sample_name = args.sample_name

    tmb_info = {
        "cov_genome": 0,
        "cov_coding": 0,
        "variants_tot": 0,
        "variants_coding": 0,
        "TMB": 0,
        "TMB_clonal": 0,
        "TMB_coding": 0,
        "TMB_clonal_coding": 0,
    }
    csin_info = {"MHCI": 0, "MHCII": 0, "combined": 0}

    tmb_info = parse_tmb(tmb, tmb_info, "all")
    tmb_info = parse_tmb(tmb_coding, tmb_info, "coding")
    csin_info = parse_csin(csin, csin_info)

    write_output(out, tmb_info, csin_info, sample_name)
    out.close()
