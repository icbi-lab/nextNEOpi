#!/usr/bin/env python

"""
This script processes protein sequences from a FASTA file, generating overlapping peptides of specified lengths.
It then outputs these peptides, along with their surrounding context sequences, to a TSV file and a FASTA file.

Requirements:
    * Python >= 3.6.2

Copyright (c) 2025 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
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
    """
    Generates sliding windows (peptides) of a specified size from a sequence and their corresponding context sequences.

    Args:
        fseq (str): The input sequence (string).
        window_size (int): The size of the sliding window (peptide length). Defaults to 8.

    Yields:
        tuple: A tuple containing:
            - win (str): The peptide sequence (string).
            - context_seq (str): The context sequence (string) surrounding the peptide.
    """
    for i in range(len(fseq) - window_size + 1):
        win = fseq[i : i + window_size]

        # Generate context sequence: 3 amino acids upstream and 3 downstream
        upstream_context = fseq[max(0, i - 3):i]  # Handle edge cases for upstream (start of sequence)
        downstream_context = fseq[i + window_size:min(len(fseq), i + window_size + 3)] # Handle edge cases for downstream (end of sequence)

        context_seq = (
            upstream_context.ljust(3, "-") +  # Pad with '-' if shorter than 3 with '-'
            win[:3] + # first 3 amino acids of the peptide
            win[-3:] + # last 3 amino acids of the peptide
            downstream_context.rjust(3,"-")  # Pad with '-' if shorter than 3 with '-'
        )
        yield win, context_seq


def chop_seqs(fasta_in, peptide_tsv_out, pep_len, fasta_out):
    """
    Chops up sequences from a FASTA file into peptides of specified lengths and writes them to a TSV file and a FASTA file.

    This function reads a FASTA file, identifies wild-type (WT) and mutant (MT) sequences,
    generates overlapping peptides of specified lengths from each sequence, and then writes
    the mutant peptides (and their context) to a TSV file and a FASTA file.

    Args:
        fasta_in (file): An open file handle for the input FASTA file.
        peptide_tsv_out (file): An open file handle for the output TSV file (peptide and context).
        pep_len (list): A list of peptide lengths to generate (e.g., [8, 9, 10, 11]).
        fasta_out (file): An open file handle for the output FASTA file (peptides).
    """
    # Group the FASTA file into headers and sequences
    faiter = (x[1] for x in groupby(fasta_in, lambda line: str(line)[0] == ">"))
    wt_seq = {}  # Dictionary to store wild-type sequences (peptide_name: peptide)
    mt_seq = {}  # Dictionary to store mutant sequences (peptide_name: peptide)
    mt_seq_context = {} # Dictionary to store mutant sequences context (peptide_name: context)

    for header in faiter:
        headerStr = str(header.__next__())
        long_name = headerStr.strip().replace(">", "")
        name = long_name.split()[0]  # Extract the sequence name (e.g., WT.gene1)
        pep_type = "wt" if name.startswith("WT") else "mt"  # Determine if it's a WT or MT sequence
        common_name = re.sub(r"^WT\.|^MT\.", "", name)  # Remove WT. or MT. prefix
        seq = "".join(str(s).strip() for s in faiter.__next__())  # Concatenate sequence lines

        i = 0
        for l in pep_len: # iterate over the peptide lengths
            for pep, context in window(seq, l): # iterate over the peptides and their context
                pep_name = common_name + ":" + str(l) + "_" + str(i) # create a unique name for the peptide
                if pep_type == "wt":
                    wt_seq[pep_name] = pep
                else:
                    mt_seq[pep_name] = pep
                    mt_seq_context[pep_name] = context
                i += 1

    wt_pep_names = set(wt_seq.keys())
    for pep_name in mt_seq:
        if pep_name in wt_pep_names:
            if mt_seq[pep_name] == wt_seq[pep_name]:
                continue # skip if the peptide is the same in wt and mt
            else:
                # peptide_tsv_out.write(">" + pep_name + "\n" + mt_seq[pep_name] + "\n")
                peptide_tsv_out.write(mt_seq[pep_name] + "\t" + mt_seq_context[pep_name] + "\n") # write the peptide and its context to the tsv file
                fasta_out.write(">" + pep_name + "\n" + mt_seq[pep_name] + "\n") # write the peptide to the fasta file


def chop_seqs_pVACsplice(fasta_in, peptide_tsv_out, pep_len, fasta_out):
    """
    Chops up sequences from a FASTA file into peptides of specified lengths and writes them to a TSV file and a FASTA file.
    This version does not discriminate between WT and MT sequences. All peptides are processed.

    Args:
        fasta_in (file): An open file handle for the input FASTA file.
        peptide_tsv_out (file): An open file handle for the output TSV file (peptide and context).
        pep_len (list): A list of peptide lengths to generate (e.g., [8, 9, 10, 11]).
        fasta_out (file): An open file handle for the output FASTA file (peptides).
    """
    # Group the FASTA file into headers and sequences
    faiter = (x[1] for x in groupby(fasta_in, lambda line: str(line)[0] == ">"))
    all_seq = {}  # Dictionary to store all sequences (peptide_name: peptide)
    all_seq_context = {} # Dictionary to store all sequences context (peptide_name: context)

    for header in faiter:
        headerStr = str(header.__next__())
        long_name = headerStr.strip().replace(">", "")
        common_name = long_name.split()[0]  # Extract the sequence name (e.g., LDHB.ENST00000350669.JUNC00000864.chr12:21654541-21654542.A)
        seq = "".join(str(s).strip() for s in faiter.__next__())  # Concatenate sequence lines

        i = 0
        for l in pep_len: # iterate over the peptide lengths
            for pep, context in window(seq, l): # iterate over the peptides and their context
                pep_name = common_name + ":" + str(l) + "_" + str(i) # create a unique name for the peptide
                all_seq[pep_name] = pep
                all_seq_context[pep_name] = context
                i += 1

    for pep_name in all_seq:
        peptide_tsv_out.write(all_seq[pep_name] + "\t" + all_seq_context[pep_name] + "\n") # write the peptide and its context to the tsv file
        fasta_out.write(">" + pep_name + "\n" + all_seq[pep_name] + "\n") # write the peptide to the fasta file


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
        "--peptide_tsv_out", required=True, type=_file_write, help="TSV file with chopped peptides ready for mixMHC2pred"
    )
    parser.add_argument(
        "--fasta_out", required=True, type=_file_write, help="Fasta file with chopped peptides"
    )
    parser.add_argument(
        "--pep_len",
        required=False,
        nargs="+",
        type=int,
        default=[8, 9, 10, 11],
        help="peptide length(s) to produce, default [8,9,10,11]",
    )
    parser.add_argument(
        "--pVACsplice",
        required=False,
        action="store_true",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    if args.pVACsplice:
        chop_seqs_pVACsplice(args.fasta_in, args.peptide_tsv_out, args.pep_len, args.fasta_out)
    else:
        chop_seqs(args.fasta_in, args.peptide_tsv_out, args.pep_len, args.fasta_out)
