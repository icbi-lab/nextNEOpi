#!/usr/bin/env python
import argparse
import sys
import csv
import subprocess
import pandas as pd

if __name__ == "__main__":

    def _file_write(fname):
        """Returns an open file handle if the given filename exists."""
        return open(fname, 'w')

    parser = argparse.ArgumentParser(description="Tanslate HLADH alleles to MixMHC2Pred alleles")
    parser.add_argument("--HLAHD_file", required=True, type=str, help="HLAHD typing file")
    parser.add_argument("--Allele_file", required=True, type=str, help="supported allele list")
    parser.add_argument("--translated_file", required=True, type=_file_write, help="translated allele list")

    args = parser.parse_args()

    #Correspondence table(dictionary) of supported alleles and corresponding nomenclatures
    with open(args.Allele_file) as file0:
        f0=csv.reader(file0,delimiter="\t")
        next(f0)
        next(f0)
        allele_dct={}
        for row0 in f0:
            allele_dct[row0[1]] =  row0[0]

    supported_alleles=list(allele_dct.keys())

    #Choose the available alleles from the result of HLA-HD
    error_colums=[]
    try:
        df=pd.read_csv(args.HLAHD_file,  sep="\t", index_col=0, header=None)
    except:
        df=pd.read_csv(args.HLAHD_file,  sep="\t", index_col=0, header=None,  error_bad_lines=False)
        with open(args.HLAHD_file, "r") as file1:
            f1=csv.reader(file1, delimiter="\t")
            for row1 in f1:
                if(len(row1) != 3):
                    error_colums.extend(row1[1:])

    allele_list_tmp=list(df[1])
    allele_list_tmp.extend(list(df[2]))
    allele_list_tmp.extend(error_colums)

    allele_list = [i for i in allele_list_tmp if (i != "-") or (i != "Not typed") ]

    for key in allele_list:
        if(key.count(":") == 2):
            if(key[:-3] in supported_alleles):
                args.translated_file.write(allele_dct[key[:-3]] + "\n")
        else:
            if(key.count(":") == 1):
                if(key in supported_alleles):
                    args.translated_file.write(allele_dct[key] + "\n")

