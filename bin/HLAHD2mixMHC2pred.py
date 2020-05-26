#!/usr/bin/python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = (
    "0",
    "2",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""

import argparse, sys

def get_hlas(in_file, hlas=[]):
    hl =[]
    hla =[]
    with open(in_file, "r") as ifile:
        for line in ifile:
            hl = line.split("\t")
            for i in range(1, len(hl)):
                hla.append(hl[i].replace("\n", ""))
        for h in hla:
            if h == "-" or h == "Not typed":
                pass
            else:
                hlas.append(h)
    return hlas

def filter_class_I(hlas=[]):
    class_I_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-V"]
    class_I_types = ["A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V"]
    class_II = []
    for hla in list(hlas):
        hla_gene = hla.split("*")[0]
        if hla_gene in class_I_genes:
            continue
        elif hla_gene not in class_I_genes and hla_gene in class_I_types:
            continue
        else:
            class_II.append(hla.replace("HLA-", ""))

    return class_II

def out_mhcii(out_file1, out_file2, hlas, supported, secondary):
    sanitized_hlas = []
    hla_ii = []
    unsupported = []
    secondary_ii = []
    # sanitize HLA list
    for hla in hlas:
        digits = hla.split(":")
        if len(digits) == 3:
            digits.pop()
            hla = (":").join(digits)
            hla = hla.replace("*", "_").replace(":", "_")
        else:
            hla = hla.replace("*", "_").replace(":", "_")
        if hla.endswith("\n"):
            hla = hla.replace("\n", "")
        sanitized_hlas.append(hla)
    # Get secondary HLAs from mixMHCpred supported HLAs
    with open(secondary) as sec:
        for row in sec:
            secondary += row.replace("\n", "")
    # Check for secondary HLAs and remove them from initial list
    for hla in sanitized_hlas:
        if hla in secondary:
            secondary_ii.append(hla)
            sanitized_hlas.remove(hla)
    secondary_ii = ', '.join(secondary_ii)

    if "DPA1_01_03" in secondary_ii and "DPA1_02_01" in secondary_ii or "DPB1_03_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_03_01")
    elif "DPA1_01_03" in secondary_ii and "DPA1_02_01" in secondary_ii or "DPB1_11_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_11_01")
    elif "DPA1_01_03" in secondary_ii and "DPA1_02_01" in secondary_ii or "DPB1_17_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_17_01")
    elif "DPA1_01_03" in secondary_ii and "DPA1_02_01" in secondary_ii or "DPB1_05_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPA1_02_02__DPB1_05_01")
    elif "DPA1_01_03" in secondary_ii or "DPB1_02_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPB1_02_01")
    elif "DPA1_01_03" in secondary_ii or "DPB1_03_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPB1_03_01")
    elif "DPA1_01_03" in secondary_ii or "DPB1_04_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPB1_04_01")
    elif "DPA1_01_03" in secondary_ii or "DPB1_06_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPB1_06_01")
    elif "DPA1_01_03" in secondary_ii or "DPB1_105_01" in secondary_ii and "DPB1_126_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPB1_105_01__DPB1_126_01")
    elif "DPA1_01_03" in secondary_ii or "DPB1_20_01" in secondary_ii:
        hla_ii.append("DPA1_01_03__DPB1_20_01")
    if "DQA1_02_01" in secondary_ii or "DQB1_02_02" in secondary_ii:
        hla_ii.append("DQA1_02_01__DQB1_02_02")
    if "DQA1_03_01" in secondary_ii or "DQB1_03_01" in secondary_ii:
        hla_ii.append("DQA1_03_01__DQB1_03_01")
    if "DQA1_03_03" in secondary_ii or "DQB1_03_01" in secondary_ii:
        hla_ii.append("DQA1_03_03__DQB1_03_01")
    if "DQA1_05_05" in secondary_ii or "DQB1_03_01" in secondary_ii:
        hla_ii.append("DQA1_05_05__DQB1_03_01")

    with open(supported, "r") as sup:
        for row in sup:
            supported += row + ", "
        for hla in sanitized_hlas:
            if hla in supported:
                hla_ii.append(hla)
            else:
                if "Not typed" in hla or "-" in hla:
                    pass
                else:
                    unsupported.append(hla.replace("_", "*", 1).replace("_", ":"))

    hla_ii = '\n'.join(hla_ii) + "\n"
    unsupported = '\n'.join(unsupported) + "\n"

    with open(out_file1, "w") as ofile1:
        ofile1.write(hla_ii)

    with open(out_file2, "w") as ofile2:
        ofile2.write("The following HLAs are not supported by mixMHC2pred:\n")
        ofile2.write(unsupported)

if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Parse HLA-HD output and format it for mixMHC2pred")

    parser.add_argument(
        "--hlahd_list", 
        required=False, 
        help="HLAs list from HLAHD"
    )
    parser.add_argument(
        "--output_dir",
        required=True, 
        help="Path to the output directory"
    )
    parser.add_argument(
        "--sample_name", 
        required=True, 
        help="Sample name"
    )
    parser.add_argument(
        "--supported_list", 
        required=True, 
        help="Supported list"
    )
    parser.add_argument(
        "--secondary_list", 
        required=True, 
        help="Secondary list"
    )
    args = parser.parse_args()
    # Parse arguments
    in_file = args.hlahd_list
    supported = args.supported_list
    secondary = args.secondary_list
    out_file_ii = args.output_dir + args.sample_name + "_mixMHC2pred.txt"
    unsup_file = args.output_dir + args.sample_name + "_unsupported.txt"


    hlas = []
    hla_ii = []
    get_hlas(in_file, hlas)
    # print(hlas)
    hla_ii = filter_class_I(hlas)
    out_mhcii(out_file_ii, unsup_file, hla_ii, supported, secondary)
