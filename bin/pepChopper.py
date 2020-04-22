#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""


import argparse
import os
import sys
import re
from itertools import groupby


def window(fseq, window_size=8):
    for i in range(len(fseq) - window_size + 1):
        yield fseq[i : i + window_size]


def chop_seqs(fasta_in, fasta_out, pep_len):
    faiter = (x[1] for x in groupby(fasta_in, lambda line: str(line)[0] == ">"))
    wt_seq = {}
    mt_seq = {}

    for header in faiter:
        headerStr = str(header.__next__())
        long_name = headerStr.strip().replace(">", "")
        name = long_name.split()[0]
        pep_type = "wt" if name.startswith("WT") else "mt"
        common_name = re.sub(r"^WT\.|^MT\.", "", name)
        seq = "".join(str(s).strip() for s in faiter.__next__())
        i = 0
        for l in pep_len:
            for pep in window(seq, l):
                pep_name = common_name + ":" + str(l) + "_" + str(i)
                if pep_type == "wt":
                    wt_seq[pep_name] = pep
                else:
                    mt_seq[pep_name] = pep
                i += 1

    wt_pep_names = set(wt_seq.keys())
    for pep_name in mt_seq:
        if pep_name in wt_pep_names:
            if mt_seq[pep_name] == wt_seq[pep_name]:
                continue
            else:
                fasta_out.write(">" + pep_name + "\n" + mt_seq[pep_name] + "\n")


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Chopup peptide sequences")

    def _file_read(fname, mode="rt"):
        """Returns an open file handle if the given filename exists."""
        if not os.path.exists(fname):
            parser.error("File '{0}' not found.".format(fname))
        _, file_extension = os.path.splitext(fname)
        if file_extension == ".gz":
            mode = "rb"
        return open(fname, mode)

    def _file_write(fname):
        """Returns an open file handle if the given filename exists."""
        return open(fname, "w")

    parser.add_argument(
        "--fasta_in", required=True, type=_file_read, help="FASTA file produced by pVACseq generate_protein_fasta"
    )
    parser.add_argument(
        "--fasta_out", required=True, type=_file_write, help="FASTA file with chopped peptides produced by this tool"
    )
    parser.add_argument(
        "--pep_len",
        required=False,
        nargs="+",
        type=int,
        default=[8, 9, 10, 11],
        help="peptide length(s) to produce, default [8,9,10,11]",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    chop_seqs(args.fasta_in, args.fasta_out, args.pep_len)
