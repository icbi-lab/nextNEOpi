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
    parser = argparse.ArgumentParser(
        description="Parse Blast result and add matches to neoantigen table"
    )
    parser.add_argument(
        "--epitope_file",
        required=True,
        type=str,
        help="TSV file with neoepitopes from pVACseq or neoFUSE",
    )
    parser.add_argument(
        "--blast_result",
        required=True,
        type=str,
        help='Output tsv from Blast run obtained using the following formating option: -outfmt "6 qseqid sseqid qlen length nident qseq"',
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
    blast_result = args.blast_result
    epitope_caller = args.epitope_caller

    epitopes = pd.read_csv(epitope_file, sep="\t")
    peptide_columnn = {"pVACseq": "MT Epitope Seq", "NeoFuse": "Fusion_Peptide"}

    blast_hits = pd.read_csv(
        blast_result,
        sep="\t",
        header=None,
        names=[
            "peptide_id",
            "ref_match_protein_id",
            "query_len",
            "align_len",
            "number_ident",
            "peptide_seq",
        ],
    )

    blast_hits = blast_hits.loc[
        (blast_hits.query_len == blast_hits.align_len)
        & (blast_hits.query_len == blast_hits.number_ident)
    ]

    blast_hits = blast_hits.groupby("peptide_seq").agg(
        {"ref_match_protein_id": lambda id: ";".join([i.split("|")[1] for i in id])}
    )

    result = pd.merge(
        epitopes,
        blast_hits,
        how="left",
        left_on=peptide_columnn[epitope_caller],
        right_on="peptide_seq",
    )

    output_file = (
        os.path.splitext(os.path.basename(epitope_file))[0] + "_ref_match.tsv"
    )
    result.to_csv(output_file, sep="\t", index=False, na_rep="-")
