#!/usr/bin/env python
"""
Requirements:
    * Python >= 3.6.2

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

import sys
import os
import argparse
import vcf
from itertools import groupby


def parse_fasta(fasta_in):
    fasta_records = []
    faiter = (x[1] for x in groupby(fasta_in, lambda line: str(line)[0] == ">"))

    for header in faiter:
        headerStr = str(header.__next__())
        # >1.MEGF6.ENST00000356575.9.missense.457P/L:17_15
        name_fields = headerStr.strip().replace(">", "").split(".")
        faiter.__next__()

        gene_symbol = ".".join(name_fields[1:-4])
        enst = name_fields[-4] + "." + name_fields[-3]
        variant_type = name_fields[-2]
        variant = name_fields[-1].split(":")[0]

        fasta_records.append(
            {"gene_symbol": gene_symbol, "enst": enst, "variant_type": variant_type, "variant": variant}
        )

    return fasta_records


def parse_vcf(vep_vcf, sample_name, normal_name):
    vcf_record = {}
    vcf_map = {}

    vcf_reader = vcf.Reader(vep_vcf)

    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS
        ref = "".join(map(str, record.REF))
        alt = "".join(map(str, record.ALT))
        af = record.genotype(sample_name)["AF"] if "AF" in str(record.genotype(sample_name)) else "NA"
        dp = record.genotype(sample_name)["DP"] if "DP" in str(record.genotype(sample_name)) else "NA"

        raf = record.genotype(sample_name)["RAF"] if "RAF" in str(record.genotype(sample_name)) else "NA"
        rdp = record.genotype(sample_name)["RDP"] if "RDP" in str(record.genotype(sample_name)) else "NA"

        gx = record.genotype(sample_name)["GX"] if "GX" in str(record.genotype(sample_name)) else "None"
        if str(gx).find("|") != -1:
            gx = str(gx).split("|")[1]

        normal_af = record.genotype(normal_name)["AF"]
        normal_dp = record.genotype(normal_name)["DP"]

        CSQ = record.INFO["CSQ"][0].split("|")
        variant_type = CSQ[1]
        gene_symbol = CSQ[3]
        ensg = CSQ[4]
        enst = CSQ[6]
        variant = CSQ[14] + CSQ[15]

        gene_symbol = gene_symbol if gene_symbol else ensg

        vcf_record = {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "af": af,
            "dp": dp,
            "raf": raf,
            "rdp": rdp,
            "gx": gx,
            "normal_af": normal_af,
            "normal_dp": normal_dp,
            "gene_symbol": gene_symbol,
            "enst": enst,
            "variant_type": variant_type,
            "variant": variant,
        }


        if "frameshift_variant" in variant_type.split("&"):
            key = gene_symbol + "_" + enst + "_" + CSQ[14] + ref + "/" + alt
        else:
            key = gene_symbol + "_" + enst + "_" + variant
        vcf_map[key] = vcf_record

    return vcf_map


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

    parser.add_argument("--sample_name", required=True, type=str, help="Tumor sample name")
    parser.add_argument("--normal_name", required=False, default="", type=str, help="Normal sample name")
    parser.add_argument("--vep_vcf", required=True, type=_file_read, help="VEP annotated VCF")
    parser.add_argument("--pep_fasta", required=True, type=_file_read, help="FASTA file with chopped peptides")
    parser.add_argument("--mixMHC2pred_result", required=True, type=_file_read, help="mixMHC2pred result file")
    parser.add_argument("--out", required=True, type=_file_write, help="Parsed mixMHC2pred result")

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    vcf_map = parse_vcf(args.vep_vcf, args.sample_name, args.normal_name)
    peptide_records = parse_fasta(args.pep_fasta)

    mixMHC2pred_result = args.mixMHC2pred_result
    out_file = args.out

    i = 0
    for line in mixMHC2pred_result:
        if line.startswith("#"):
            continue
        if line.startswith("Peptide"):
            out_file.write(
                "Chromosome"
                + "\t"
                + "Positions"
                + "\t"
                + "Ref"
                + "\t"
                + "Alt"
                + "\t"
                + "AF"
                + "\t"
                + "DP"
                + "\t"
                + "RAF"
                + "\t"
                + "RDP"
                + "\t"
                + "GX"
                + "\t"
                + "Normal_AF"
                + "\t"
                + "Normal_DP"
                + "\t"
                + "Symbol"
                + "\t"
                + "Transcript"
                + "\t"
                + "Variant_Type"
                + "\t"
                + "Variant"
                + "\t"
                + line
            )
            continue
        key = peptide_records[i]["gene_symbol"] + "_" + peptide_records[i]["enst"] + "_" + peptide_records[i]["variant"]
        out_file.write("\t".join(map(str, vcf_map[key].values())) + "\t" + line)
        i += 1
