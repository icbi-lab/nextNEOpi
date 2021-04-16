#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2
    * Pandas
    * NumPy

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
import math
import os, sys
import pandas as pd
import numpy as np

if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Add CCF to neoepitopes")

    parser.add_argument("--neoepitopes", required=True, help="neoepitope file")
    parser.add_argument("--ccf", required=True, help="CCF file")
    parser.add_argument("--outfile", required=True, help="Output file")

    args = parser.parse_args()

    neoepitope_reader = pd.read_csv(args.neoepitopes, sep="\t")
    neoepitope_df = pd.DataFrame(neoepitope_reader)

    ccf_reader = pd.read_csv(args.ccf, sep="\t")
    ccf_df = pd.DataFrame(ccf_reader)

    # make positions 0 based
    ccf_df["pos"] -= 1

    neoepitope_df = pd.merge_ordered(
        neoepitope_df,
        ccf_df[["Chr", "pos", "CCF", "CCF.05", "CCF.95", "pSubclonal", "pClonal"]],
        how="left",
        left_on=["Chromosome", "Start"],
        right_on=["Chr", "pos"],
        fill_method="ffill",
    ).drop(["Chr", "pos"], axis=1)

    neoepitope_df.to_csv(args.outfile, sep="\t", index=False)
