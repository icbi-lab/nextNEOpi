#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2
    * Pysam

Copyright (c) 2021 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>


"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""


import pysam
import argparse
import sys
import os
import numpy as np
from multiprocessing import Pool


# Global
coding_variants = [
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
]


def is_coding(variant_type):
    v_types = variant_type.split("&")
    for v_type in v_types:
        if v_type in coding_variants:
            return True

    return False


def get_coverage(tumor_bam, normal_bam, regions, mode, min_q, min_c, fasta_file):
    coverage = 0
    truncate = True if mode == "truncate" else False

    bam_t_fh = pysam.pysam.AlignmentFile(tumor_bam, "rb", threads = 2)
    bam_n_fh = pysam.pysam.AlignmentFile(normal_bam, "rb", threads = 2)

    fasta = None
    if fasta_file:
        try:
            fasta = pysam.FastaFile(fasta_file)
        except ValueError:
            print("Please use samtools to create a faidx for: " + fasta_file)
            raise
        except IOError:
            print("File not found: " + fasta_file)
            raise


    for region in regions:
        # print(region)
        contig, start, end = region

        region_size = end - start + 1
        coverage_t = np.full(region_size, False)
        coverage_n = np.full(region_size, False)

        for pileupcolumn_t in bam_t_fh.pileup(contig=contig, start=start, stop=end, compute_baq = False, min_base_quality=min_q, truncate=truncate, stepper="samtools", fastafile=fasta):
            n = pileupcolumn_t.get_num_aligned()
            if n >= min_c:
                coverage_t[pileupcolumn_t.pos - start + 1] = True

        for pileupcolumn_n in bam_n_fh.pileup(contig=contig, start=start, stop=end, compute_baq = False, min_base_quality=min_q, truncate=truncate, stepper="samtools", fastafile=fasta):
            n = pileupcolumn_n.get_num_aligned()
            if n >= min_c:
                coverage_n[pileupcolumn_n.pos - start + 1] = True

        coverage += np.sum(np.multiply(coverage_t, coverage_n))

    return coverage


def get_coding_regions_from_bed(bed_file, chrom):
    coding_regions = list()

    if not bed_file.endswith(".gz") and not os.path.exists(bed_file + ".gz"):
        try:
            bgz_bed = pysam.tabix_index(bed_file, preset="bed", force=True, keep_original=True)
        except IOError:
            print("File not found: " + bed_file)
            raise

    elif bed_file.endswith(".gz") and not os.path.exists(bed_file + ".gz.tbi"):
        try:
            bgz_bed = pysam.tabix_index(bed_file, preset="bed", force=True, keep_original=True)
        except IOError:
            print("File not found: " + bed_file)
            raise
    else:
        bgz_bed = bed_file + ".gz"

    bed = pysam.TabixFile(bgz_bed, parser=pysam.asBed())
    for row in bed.fetch(chrom):
            coding_regions.append([row.contig, row.start, row.end])

    return coding_regions


def parse_vcf_info(rec):
    INFO = {}
    info_recs = rec.info.split(";")
    for info_rec in info_recs:
        if '=' in info_rec:
            k, v = info_rec.split("=")
            INFO[k] = v
        else:
            INFO[info_rec] = ""
    return INFO


def get_variants_from_vcf(vcf_file, var_type):
    variants = 0

    if not vcf_file.endswith(".gz") and not os.path.exists(vcf_file + ".gz"):
        try:
            bgz_vcf = pysam.tabix_index(vcf_file, preset="vcf", force=True, keep_original=True)
        except IOError:
            print("File not found: " + vcf_file)
            raise

    elif vcf_file.endswith(".gz") and not os.path.exists(vcf_file + ".gz.tbi"):
        try:
            bgz_vcf = pysam.tabix_index(vcf_file, preset="vcf", force=True, keep_original=True)
        except IOError:
            print("File not found: " + vcf_file)
            raise
    else:
        bgz_vcf = vcf_file + ".gz"

    vcf = pysam.TabixFile(bgz_vcf, parser=pysam.asVCF())

    for rec in vcf.fetch():
        if var_type == "coding":
            INFO = parse_vcf_info(rec)
            CSQ = INFO["CSQ"].split("|")
            variant_type = CSQ[1]
            if is_coding(variant_type):
                variants += 1
        else:
            variants += 1

    return variants



def _file_write(fname):
    """Returns an open file handle if the given filename exists."""
    return open(fname, "w")


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        description="Calculate the tumor mutational burden or load"
    )
    parser.add_argument(
        "--normal_bam",
        required=True,
        type=str,
        help="BAM file from normal tissue",
    )
    parser.add_argument(
        "--tumor_bam",
        required=True,
        type=str,
        help="BAM file from tumor tissue",
    )
    parser.add_argument(
        "--vcf",
        required=True,
        type=str,
        help="VCF file with variants from tumor",
    )
    parser.add_argument(
        "--bed",
        required=False,
        type=str,
        default="",
        help="BED file with merged exons.",
    )
    parser.add_argument(
        "--fasta",
        required=False,
        type=str,
        help="Fasta file of the genome. Same as used to generate the BAM files",
    )
    parser.add_argument(
        "--output_file",
        required=True,
        type=_file_write,
        help="Output file storing the results",
    )
    parser.add_argument(
        "--min_BQ",
        required=False,
        type=int,
        default=20,
        help="Minimum basecall quality to consider",
    )
    parser.add_argument(
        "--min_coverage",
        required=False,
        type=int,
        default=10,
        help="Minimum read coverage of a variant positon",
    )
    parser.add_argument(
        "--variant_type",
        required=False,
        type=str,
        default="all",
        choices=["all", "coding"],
        help="Variant types to consider",
    )
    parser.add_argument(
        "--cpus",
        required=False,
        type=int,
        default=1,
        help="Number of CPUs to use for parallel processing.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__
    )

    args = parser.parse_args()

    tumor_bam = args.tumor_bam
    normal_bam = args.normal_bam
    min_q = args.min_BQ
    min_c = args.min_coverage
    vcf_file = args.vcf
    fasta_file = args.fasta
    bed_file = args.bed
    var_type = args.variant_type
    output_file = args.output_file
    cpus = args.cpus

    variants = get_variants_from_vcf(vcf_file, var_type)

    bam = pysam.AlignmentFile(tumor_bam, "rb", threads = 2)

    chr_start = 1
    mode = "all"

    pool = Pool(processes=cpus)
    pool_result = {}
    coding_regions = {}

    for chrom in bam.references:
        # skip short contigs (chr 21 is shortest ~45 Mbases)
        if bam.get_reference_length(chrom) > 40000000:
            pool_result[chrom] = None
            if bed_file is not "":
                coding_regions[chrom] = get_coding_regions_from_bed(bed_file, chrom)
                mode = "truncate"


    for chrom in pool_result.keys():
        chrom_len = bam.get_reference_length(chrom)

        if len(coding_regions) == 0:
            regions = [[chrom, chr_start, chr_start + chrom_len - 1]]
        else:
            regions = coding_regions[chrom]

        pool_result[chrom] = pool.apply_async(get_coverage, (tumor_bam, normal_bam, regions, mode, min_q, min_c, fasta_file))

    pool.close()
    pool.join()

    coverage = np.sum([pool_result[chrom].get() for chrom in pool_result.keys()])
    del pool_result

    genome_part = "entire genome" if bed_file == "" else "exons"
    output_file.write("Coverage over " + genome_part + ":\t" + str(coverage) + "\n")
    output_file.write("Min depth in T and N:\t" + str(min_c) + "\n")
    output_file.write("Min base call quality considerd:\t" + str(min_q) + "\n")
    output_file.write("Variants of type " + var_type + ":\t" + str(variants) + "\n")
    output_file.write("Mutational load (variants/Mbps):\t" + str(round(variants*10**6/coverage, 3)) + "\n")
