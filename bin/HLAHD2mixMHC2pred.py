#!/usr/bin/env python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = (
    "0",
    "3",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""

import argparse, sys


def get_hlas(in_file, hlas=[]):
    hl = []
    hla = []
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
            if "Not" in hla or "typed" in hla or hla == "":
                pass
            else:
                class_II.append(hla.replace("HLA-", ""))

    return class_II


def out_mhcii(out_file1, out_file2, out_file3, hlas, supported, model_list):
    sanitized_hlas = []
    # DRB = []
    DPA = []
    DPB = []
    DPB = []
    DQA = []
    DQB = []
    supp = []
    scnd = []
    model = []
    hla_ii = []
    hla_conf = []
    unsupported = []
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

    # Get supported HLAs
    with open(supported) as sup:
        for row in sup:
            supp.append(row.replace("\n", ""))
    # Get mixMHC2pred models
    with open(model_list) as mod:
        for row in mod:
            model.append(row.replace("\n", ""))

    # Get valid predicted HLA types (and process DRB alleles)
    for hla in sanitized_hlas:
        if "DRB" in hla:
            if hla in supp:
                hla_conf.append(hla + "\tH")
                hla_ii.append(hla)
            else:
                unsupported.append(hla)
        elif "DP" in hla:
            if "DPA" in hla and hla in supp:
                DPA.append(hla)
            elif "DPA" in hla and hla not in supp:
                unsupported.append(hla)
            elif "DPB" in hla and hla in supp:
                DPB.append(hla)
            elif "DPB" in hla and hla not in supp:
                unsupported.append(hla)
        elif "DQ" in hla:
            if "DQA" in hla and hla in supp:
                DQA.append(hla)
            elif "DQA" in hla and hla not in supp:
                unsupported.append(hla)
            elif "DQB" in hla and hla in supp:
                DQB.append(hla)
            elif "DQB" in hla and hla not in supp:
                unsupported.append(hla)
    # Process DQ
    ## Find all possible combos of DQ
    DQ_AB = []
    if DQA and not DQB:
        unsupported.append(DQA[0])
    elif not DQA and DQB:
        DQ_AB.append(DQB[0])
    else:
        for dqa in DQA:
            for dqb in DQB:
                DQ_AB.append(dqa + "__" + dqb)

    ## Run HLA-DQ combos against known models
    for dq in DQ_AB:
        if dq == "DQB1_02_02":
            hla_conf.append("DQA1_02_01__DQB1_02_02\tL")
            hla_ii.append("DQA1_02_01__DQB1_02_02")
        elif dq == "DQB1_03_01":
            hla_conf.append("DQA1_03_01__DQB1_03_01\tL")
            hla_ii.append("DQA1_03_01__DQB1_03_01")
        else:
            for mod in model:
                if dq == mod:
                    hla_conf.append(dq + "\tH")
                    hla_ii.append(dq)

    # Process DP
    ## Find all possible combos of DQ
    DP_AB = []
    if DPA and not DPB:
        if "DPA1_01_03" in DPA:
            if "DPA1_02_01" in DPA:
                hla_conf.append("DPA1_01_03__DPA1_02_01__DPB1_03_01\tL")
                hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_03_01")
            else:
                unsupported.append("DPA1_01_03")
    elif not DPA and DPB:
        if "DPB1_105_01" in DPB:
            if "DPB1_126_01" in DPB:
                hla_conf.append("DPA1_01_03__DPB1_105_01__DPB1_126_01\tL")
                hla_ii.append("DPA1_01_03__DPB1_105_01__DPB1_126_01")
                DPB.remove("DPB1_105_01")
                DPB.remove("DPB1_126_01")
        else:
            for dpb in DPB:
                for mod in model:
                    if dpb in mod:
                        hla_conf.append(mod + "\tL")
                        hla_ii.append(mod)

    if DPA and DPB:
        if "DPA1_01_03" in DPA and "DPA1_02_01" in DPA:
            if "DPB1_03_01" not in DPB:
                hla_conf.append("DPA1_01_03__DPA1_02_01__DPB1_03_01\tL")
                hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_03_01")
            elif "DPB1_11_01" not in DPB:
                hla_conf.append("DPA1_01_03__DPA1_02_01__DPB1_03_01\tL")
                hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_03_01")
            elif "DPB1_17_01" not in DPB:
                hla_conf.append("DPA1_01_03__DPA1_02_01__DPB1_03_01\tL")
                hla_ii.append("DPA1_01_03__DPA1_02_01__DPB1_03_01")
            else:
                DPA.append("DPA1_01_03__DPA1_02_01")
                DPA.remove("DPA1_02_01")
        elif "DPA1_01_03" in DPA and "DPA1_02_02" in DPA:
            if "DPB1_05_01" not in DPB:
                hla_conf.append("DPA1_01_03__DPA1_02_02__DPB1_05_01\tL")
                hla_ii.append("DPA1_01_03__DPA1_02_02__DPB1_05_01")
            else:
                DPA.append("DPA1_01_03__DPA1_02_02")
                DPA.remove("DPA1_02_02")
        if "DPB1_105_01" in DPB:
            if "DPB1_126_01" in DPB:
                DPB.append("DPB1_105_01__DPB1_126_01")
                DPB.remove("DPB1_105_01")
                DPB.remove("DPB1_126_01")
        for dpa in DPA:
            for dpb in DPB:
                DP_AB.append(dpa + "__" + dpb)
    for dp in DP_AB:
        if dp == "DPA1_01_03__DPB1_105_01" or dp == "DPA1_01_03__DPB1_126_01":
            hla_conf.append("DPA1_01_03__DPB1_105_01__DPB1_126_01\tH")
            hla_ii.append("DPA1_01_03__DPB1_105_01__DPB1_126_01")
        else:
            for mod in model:
                if dp == mod:
                    hla_conf.append(mod + "\tH")
                    hla_ii.append(mod)

    with open(out_file1, "w") as ofile1:
        hla_ii_no_dups = list(dict.fromkeys(hla_ii))
        ofile1.write("\n".join(hla_ii_no_dups))

    with open(out_file2, "w") as ofile2:
        ofile2.write("The following HLAs are not supported by mixMHC2pred:\n")
        unsupported_no_dups = list(dict.fromkeys(unsupported))
        ofile2.write("\n".join(unsupported_no_dups))

    with open(out_file3, "w") as ofile3:
        hla_conf_no_dups = list(dict.fromkeys(hla_conf))
        ofile3.write("\n".join(hla_conf_no_dups))


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Parse HLA-HD output and format it for mixMHC2pred")

    parser.add_argument("--hlahd_list", required=False, help="HLAs list from HLAHD")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    parser.add_argument("--supported_list", required=True, help="Supported list")
    parser.add_argument("--model_list", required=True, help="Secondary list")
    args = parser.parse_args()
    # Parse arguments
    in_file = args.hlahd_list
    supported = args.supported_list
    model_list = args.model_list
    out_file_ii = args.output_dir + args.sample_name + "_mixMHC2pred.txt"
    unsup_file = args.output_dir + args.sample_name + "_unsupported.txt"
    confidence_file = args.output_dir + args.sample_name + "_mixMHC2pred_conf.txt"

    hlas = []
    hla_ii = []
    get_hlas(in_file, hlas)
    hla_ii = filter_class_I(hlas)
    out_mhcii(out_file_ii, unsup_file, confidence_file, hla_ii, supported, model_list)
