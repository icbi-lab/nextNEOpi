#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2
    * Pandas

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


import pandas as pd
import argparse
import sys
import os
import gzip


def _file_write(fname):
    """Returns an open file handle if the given filename exists."""
    return open(fname, "w")


def _file_read(fname, mode="rt"):
    """Returns an open file handle if the given filename exists."""
    if not os.path.exists(fname):
        parser.error("File '{0}' not found.".format(fname))
    _, file_extension = os.path.splitext(fname)
    if file_extension == ".gz":
        mode = "rb"
        return gzip.open(fname, mode), mode
    else:
        return open(fname, mode), mode


if __name__ == "__main__":
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Make peptide fasta")
    parser.add_argument(
        "--epitope_file",
        required=True,
        type=str,
        help="TSV file with neoepitopes from pVACseq or neoFUSE",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=_file_write,
        help="Output fasta file",
    )
    parser.add_argument(
        "--epitope_caller",
        required=True,
        type=str,
        choices=["pVACseq", "NeoFuse"],
        help="Epitope calling tool used: [pVACseq|NeoFuse]",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args()

    epitope_file = args.epitope_file
    fasta = args.fasta
    epitope_caller = args.epitope_caller

    epitopes = pd.read_csv(epitope_file, sep="\t")

    peptide_columnn = {"pVACseq": "MT Epitope Seq", "NeoFuse": "Fusion_Peptide"}

    idx = 0
    for pep_seq in epitopes[peptide_columnn[epitope_caller]].unique():
        fasta.write(">pep_" + str(idx) + "\n" + pep_seq + "\n")
        idx += 1
