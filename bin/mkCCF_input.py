#!/usr/bin/env python

"""
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


# Global
coding_variants = [
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "coding_sequence_variant",
]


def is_number(a):
    try:
        float(a)
        bool_a = True
    except:
        bool_a = False

    return bool_a


def get_segments(seg_file):
    segments = {}

    i = 0
    old_chrom = ""

    # skip headr
    seg_file.readline()

    for line in seg_file:
        line = line.strip()
        chrom, start_pos, end_pos, n_major, n_minor = line.split("\t")
        if chrom != old_chrom:
            i = 0
            old_chrom = chrom
            segments["chr_prefix"] = True if chrom.startswith("chr") else False
            segments[chrom] = []

        if n_major == "NA" or n_minor == "NA":
            continue

        segments[chrom].append(
            {
                "start_pos": int(float(start_pos)),
                "end_pos": int(float(end_pos)),
                "n_major": int(n_major),
                "n_minor": int(n_minor),
            }
        )

    seg_file.close()

    return segments


def get_segments_seqz(seg_file):
    segments = {}

    i = 0
    old_chrom = ""

    # skip headr
    seg_file.readline()

    for line in seg_file:
        line = line.strip()
        fields = line.split("\t")
        chrom, start_pos, end_pos, n_major, n_minor = [fields[0], fields[1], fields[2], fields[10], fields[11]]
        if chrom != old_chrom:
            i = 0
            old_chrom = chrom
            segments["chr_prefix"] = True if chrom.startswith("chr") else False
            segments[chrom] = []

        if n_major == "NA" or n_minor == "NA":
            continue

        segments[chrom].append(
            {
                "start_pos": int(float(start_pos)),
                "end_pos": int(float(end_pos)),
                "n_major": int(n_major),
                "n_minor": int(n_minor),
            }
        )

    seg_file.close()

    return segments


def get_purity(purity_file):
    for line in purity_file:
        line = line.strip().split("\t")[0]

    purity_file.close()

    return line


def get_purity_seqz(purity_file):

    purity_file.readline()
    purity_file.readline()

    line = purity_file.readline()
    purity = line.strip().split("\t")[0]

    purity_file.close()

    return purity


def is_coding(variant_type):
    v_types = variant_type.split("&")
    for v_type in v_types:
        if v_type in coding_variants:
            return True

    return False


def make_ccf_calc_input(pat_id, sample_type, vcf_in, segments, purity, min_vaf, write_fh):
    sampleName = pat_id
    reader = vcf.Reader(vcf_in)
    header_fields = [
        "PatientID",
        "SampleType",
        "Chr",
        "pos",
        "normCnt",
        "mutCnt",
        "VAF",
        "ploidy",
        "cnA",
        "cnB",
        "purity",
        "coding",
    ]
    write_fh.write("\t".join(header_fields) + "\n")

    for rec in reader:
        vaf = rec.genotype(sampleName)["AF"]
        ad = rec.genotype(sampleName)["AD"]
        norm_cnt, mut_cnt = ad

        variant_type = rec.INFO["CSQ"][0].split("|")[1]
        coding = "T" if is_coding(variant_type) else "F"

        chrom_key = (
            rec.CHROM.lstrip("chr") if (rec.CHROM.startswith("chr")) and not segments["chr_prefix"] else rec.CHROM
        )

        # TODO: deal with chrM
        if chrom_key not in segments:
            continue

        # set defaults
        seg_found = False
        ploidy = 2
        cn_a = 1
        cn_b = 1

        for seg in segments[chrom_key]:
            if rec.POS >= seg["start_pos"] and rec.POS <= seg["end_pos"]:
                ploidy = seg["n_minor"] + seg["n_major"]
                cn_a = seg["n_major"]
                cn_b = seg["n_minor"]
                seg_found = True
                break

        if not seg_found:
            print("WARNING: " + rec.CHROM + ":" + str(rec.POS) + " not within any CN segment: setting cn_minor and cn_major to 1")

        rec_out = [pat_id, sample_type, rec.CHROM, rec.POS, norm_cnt, mut_cnt, vaf, ploidy, cn_a, cn_b, purity, coding]

        if vaf >= min_vaf:
            write_fh.write("\t".join(map(str, rec_out)) + "\n")

    write_fh.close()


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Create vcf with coding only variants")

    def _file_read(fname, mode="rt"):
        """Returns an open file handle if the given filename exists."""
        if not os.path.exists(fname):
            parser.error("File '{0}' not found.".format(fname))
        _, file_extension = os.path.splitext(fname)
        if file_extension == ".gz":
            mode = "rb"
        return open(fname, mode)

    def _file_write(fname, mode="w"):
        """Returns an open file handle if the given filename exists."""
        return open(fname, mode)

    parser.add_argument("--PatientID", required=True, default="", type=str, help="Patient ID")
    parser.add_argument("--sample_type", required=False, default="Tumor", type=str, help="Sample type string")
    parser.add_argument("--vcf", required=True, type=str, help="VCF file")
    parser.add_argument("--seg", required=False, default="", type=str, help="ASCAT segment/cnv file *.cnvs.txt")
    parser.add_argument(
        "--seg_sequenza", required=False, default="", type=str, help="Sequenza segment/cnv file segments.txt"
    )
    parser.add_argument(
        "--purity",
        required=False,
        type=str,
        help="Tumor purity as nmuber or ASCAT tumor purity file *.purityploidy.txt",
    )
    parser.add_argument(
        "--purity_sequenza",
        required=False,
        type=str,
        help="Tumor purity as nmuber or sequenza tumor purity file confints_CP.txt",
    )

    parser.add_argument("--min_vaf", required=False, default=0.05, type=float, help="minimum VAF to consider")
    parser.add_argument(
        "--result_table",
        required=False,
        default="CCF_input.tsv",
        type=str,
        help="Output file for CCF calculations input",
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    if args.seg != "" and args.seg_sequenza != "":
        print("ERROR: please provide a segments file from ASCAT or Sequenza")
        sys.exit(1)

    if args.purity_sequenza:
        purity = (
            args.purity_sequenza
            if is_number(args.purity_sequenza)
            else get_purity_seqz(_file_read(args.purity_sequenza))
        )
    else:
        purity = args.purity if is_number(args.purity) else get_purity(_file_read(args.purity))

    if float(purity) > 1 or float(purity) <= 0 or is_number(purity) is False:
        print(is_number(purity))
        print("Warning: purity should be 0> p <=1, got: " + str(purity) + " - setting it to 0.75")
        purity = "0.75"

    if args.seg_sequenza:
        segments = get_segments_seqz(_file_read(args.seg_sequenza))
    else:
        segments = get_segments(_file_read(args.seg))

    make_ccf_calc_input(
        args.PatientID,
        args.sample_type,
        _file_read(args.vcf),
        segments,
        purity,
        args.min_vaf,
        _file_write(args.result_table),
    )
