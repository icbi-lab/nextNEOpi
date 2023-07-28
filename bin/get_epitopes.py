#!/usr/bin/env python

"""
Parse the output files of pVACseq to get all epitopes (WT and mutated), with the mutation position

Requirements:
    * Python >= 2.7

"""
import argparse
import os
import sys
import csv

def filter_tsv(sample, input_file, output_file):
    # Open the input TSV file and create a CSV reader
    with open(input_file, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')

        # Open the output TSV file and create a CSV writer
        with open(output_file, 'w', newline='') as outfile:
            fieldnames = ['Sample_ID', 'mut_peptide', 'Reference', 'peptide_variant_position']
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')

            # Write the header
            writer.writeheader()

            # Filter rows and write selected columns to the output TSV file
            for row in reader:
                if row['WT Epitope Seq'] and row['MT Epitope Seq']:
                    writer.writerow({'Sample_ID': sample,
                                     'mut_peptide': row['MT Epitope Seq'],
                                     'Reference': row['WT Epitope Seq'],
                                     'peptide_variant_position': row['Mutation Position']})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse the output files of pVACseq for all MHC I epitopes")
    parser.add_argument("--pvacseq_out", required=True, help="output file produced by pvacseq")
    parser.add_argument("--output", required=True, help="output file produced by tis script")
    parser.add_argument("--sample_id", required=True, help="sample name")
    args = parser.parse_args()
    infile = args.pvacseq_out
    outfile = args.output
    sample = args.sample_id
    epitope_array = []

    filter_tsv(sample, infile, outfile)