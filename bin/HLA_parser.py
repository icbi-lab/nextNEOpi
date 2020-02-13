#!/usr/bin/env python

"""
Parse the output of OptiType and HLA-HD and get HLA types (I and II) at 4 digit resolution

Requirements:
    * Python >= 3.6.2

"""
import argparse
import os
import sys
import csv

def parse_hlahd(inFile, hlas=[]):
	with open(inFile) as in_file:
		csv_reader = csv.reader(in_file, delimiter='\t')
		in_file.readline()
		for line in csv_reader:
			if "HLA-" in line[1]:
				hlas.append(line[1].split(":")[0]+":"+line[1].split(":")[1])
			if "HLA-" in line[2]:
				hlas.append(line[2].split(":")[0]+":"+line[2].split(":")[1])
	return hlas

def parse_opti(inFile, hlas=[]):
	with open(inFile) as in_file:
		csv_reader = csv.reader(in_file, delimiter='\t')
		in_file.readline()
		for line in csv_reader:
			if "*" in line[1]:
				hlas.append("HLA-"+line[1])
			if "*" in line[2]:
				hlas.append("HLA-"+line[2])
			if "*" in line[3]:
				hlas.append("HLA-"+line[3])
			if "*" in line[4]:
				hlas.append("HLA-"+line[4])
			if "*" in line[5]:
				hlas.append("HLA-"+line[5])
			if "*" in line[6]:
				hlas.append("HLA-"+line[6])
	return hlas

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parse the output files of OptiType and HLA-HD for all HLA types and alleles')
	parser.add_argument('--opti_out', required=True, help='output file produced by OptiType')
	parser.add_argument('--hlahd_out', required=True, help='output file produced by HLA-HD')
	args = parser.parse_args()
	# infile = "/data/projects/2019/NeoAG/VCF-phasing/RESULTS/CRC19/09_HLA_HD/CRC19/result/CRC19_final.result.txt"
	infile = args.opti_out
	# infile2 = "/data/projects/2019/NeoAG/VCF-phasing/RESULTS/CRC19/08_OptiType/2020_02_11_15_35_04/2020_02_11_15_35_04_result.tsv"
	infile2 = args.hlahd_out
	hla_array = []
	parse_opti(infile, hla_array)
	parse_hlahd(infile2, hla_array)
	print(",".join(hla_array))