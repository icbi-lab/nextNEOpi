#!/usr/bin/env python

"""
Merge VCF files called by GATK's HaplotypeCaller & UnifiedGenotyper
Code was inspired by Wibowo Arindrarto <w.arindrarto@lumc.nl>

Requirements:
    * Python >= 3.6.2
    * PyVCF >= 0.6.8

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
import os
import sys

import vcf
from vcf.parser import _Info
from vcf.utils import walk_together


def make_hc_somatic_vars(f_primary, primary_caller_name, f_confirming, confirming_caller_names_, fout, fout_single):
    """Picks the given VARs from the priority vcf to a new
    VCF if they are confirmed by one of the others."""

    var_count = 0

    primary_reader = vcf.Reader(f_primary)

    confirming_reader = []
    for f in f_confirming:
        confirming_reader.append(vcf.Reader(f))

    primary_reader.infos["VariantCalledBy"] = _Info(
        "VariantCalledBy", ".", "String", "variant callers that called the variant", "caller ", "0.1"
    )

    # some sanity checks
    sorted_primary_contigs = sorted(primary_reader.contigs.keys())
    for cr in confirming_reader:
        if sorted(primary_reader.samples) != sorted(cr.samples):
            raise ValueError("Input VCF files must have the same sample column " "headers.")
        if sorted_primary_contigs != sorted(cr.contigs.keys()):
            raise ValueError("Input VCF files must denote the same contigs.")

    out_writer = vcf.Writer(fout, primary_reader)
    out_writer_single = vcf.Writer(fout_single, primary_reader)

    all_readers = [primary_reader]
    for r in confirming_reader:
        all_readers.append(r)

    for primary_rec, *confirming_recs in walk_together(*all_readers):
        confirmed = False
        if primary_rec is not None:
            confirmed_idx = [i for i, x in enumerate(confirming_recs) if x is not None]
            if len(confirmed_idx) > 0:
                confirming_caller_names = ",".join(confirming_caller_names_[i] for i in confirmed_idx)
                primary_rec.add_info("VariantCalledBy", primary_caller_name + "," + confirming_caller_names)
                confirmed = True
            else:
                primary_rec.add_info("VariantCalledBy", primary_caller_name)
                out_writer_single.write_record(primary_rec)

            if confirmed:
                out_writer.write_record(primary_rec)
                var_count += 1

    return var_count


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        description="Create vcf that contains only vars confirmed by at least one other caller"
    )

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
        "--primary",
        required=True,
        type=_file_read,
        help="VCF file from which variants are kept when confirmed by any other vcf in --confirming",
    )
    parser.add_argument(
        "--primary_name", required=True, type=str, help="Name of the primary variant caller",
    )
    parser.add_argument(
        "--confirming",
        required=True,
        nargs="*",
        type=_file_read,
        help="VCF files used to confirm variants in --primary",
    )
    parser.add_argument(
        "--confirming_names", required=True, type=str, nargs="*", help="Names of the confirming variant callers",
    )
    parser.add_argument(
        "--out_vcf", required=True, type=_file_write, help="VCF file for confirmed vars produced by " "this tool"
    )
    parser.add_argument(
        "--out_single_vcf",
        required=True,
        type=_file_write,
        help="VCF file for unconfirmed vars produced by " "this tool",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    var_count = make_hc_somatic_vars(
        args.primary, args.primary_name, args.confirming, args.confirming_names, args.out_vcf, args.out_single_vcf
    )

    print("Total confirmed:\t" + str(var_count))
