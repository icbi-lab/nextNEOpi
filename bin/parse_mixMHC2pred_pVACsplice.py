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

import sys
import os
import argparse
import csv
from itertools import groupby


def parse_fasta(fasta_in):
    fasta_records = []
    faiter = (x[1] for x in groupby(fasta_in, lambda line: str(line)[0] == ">"))

    for header in faiter:
        headerStr = str(header.__next__())
        # >LDHB.ENST00000350669.JUNC00000864.chr12:21654541-21654542.A:15_0
        name_fields = headerStr.strip().replace(">", "").split(".")
        faiter.__next__()

        gene_symbol = name_fields[0]
        enst = name_fields[1]
        variant = name_fields[2]
        variant_info = name_fields[3]

        fasta_records.append(
            {"gene_symbol": gene_symbol, "enst": enst, "variant": variant, "variant_info": variant_info}
        )

    return fasta_records


def create_regtools_map(regtools_tsv):
    """
    Reads a TSV file with specific columns and returns a map with keys
    constructed from "name" and "variant_info".

    Args:
        tsv_filepath (str): The path to the TSV file.

    Returns:
        dict: A dictionary where keys are strings formed by concatenating
              the "name" and "variant_info" columns with an underscore,
              and values are the corresponding rows as lists.
    """
    regtools_map = {}

    reader = csv.reader(regtools_tsv, delimiter='\t')
    l = 0

    for row in reader:
        if l == 0:
            header = row
            l += 1
            continue
        if len(row) >= 18:
            name = row[3]
            variant_info = row[17]
            key = f"{name}_{variant_info}"
            regtools_map[key] = row
        else:
            print(f"Warning: Skipping row with fewer than 18 columns: {row}")
        l += 1
    return header, regtools_map



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

    parser.add_argument("--regtools_tsv", required=True, type=_file_read, help="regtools TSV file")
    parser.add_argument("--pep_fasta", required=True, type=_file_read, help="FASTA file with chopped peptides")
    parser.add_argument("--mixMHC2pred_result", required=True, type=_file_read, help="mixMHC2pred result file")
    parser.add_argument("--out", required=True, type=_file_write, help="Parsed mixMHC2pred result")

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    header, regtools_map = create_regtools_map(args.regtools_tsv)
    peptide_records = parse_fasta(args.pep_fasta)

    mixMHC2pred_result = args.mixMHC2pred_result
    out_file = args.out

    i = 0
    for line in mixMHC2pred_result:
        if line.startswith("#"):
            continue
        if line.startswith("Peptide"):
            out_file.write(
                "\t".join(header)
                + "\t"
                + line
            )
            continue
        key = peptide_records[i]["variant"] + "_" + peptide_records[i]["variant_info"]
        out_file.write("\t".join(map(str, regtools_map[key])) + "\t" + line)
        i += 1
