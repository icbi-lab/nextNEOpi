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

def parse_mhcI(inFile, epitopes=[]):
	with open(inFile) as in_file:
		csv_reader = csv.reader(in_file, delimiter='\t')
		in_file.readline()
		for line in csv_reader:
			if line[19] == "NA" or line[18] == "NA":
				pass
			else:
			# 	print("%s\t%s\t%s" % (line[18], line[19], line[17]))
				epitopes.append("%s\t%s\t%s" % (line[18], line[19], line[17]))
	return epitopes

def write_output(outFile, sample_id, epitopes=[]):
	with open(outFile, "w") as out_file:
		out_file.write("Sample_ID\tmut_peptide\tReference\tpeptide_variant_position\n")
		for epitope in epitopes:
			out_file.write("%s\t%s\n" % (sample_id, epitope))
	return outFile

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parse the output files of pVACseq for all MHC I epitopes')
	parser.add_argument('--pvacseq_out', required=True, help='output file produced by pvacseq')
	parser.add_argument('--output', required=True, help='output file produced by tis script')
	parser.add_argument('--sample_id', required=True, help='sample name')
	args = parser.parse_args()
	infile = args.pvacseq_out
	outfile = args.output
	sample = args.sample_id
	epitope_array = []

	parse_mhcI(infile, epitope_array)
	write_output(outfile, sample, epitope_array)