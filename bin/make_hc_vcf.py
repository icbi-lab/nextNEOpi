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
from vcf.parser import _Info, _Format
from vcf.utils import walk_together


def make_hc_somatic_vars(f_primary, primary_caller_name, f_confirming, confirming_caller_names_, fout, fout_single):
    """Picks the given VARs from the priority vcf to a new
    VCF if they are confirmed by one of the others."""

    var_count = 0
    normal_sample_name = args.normal_sample_name

    primary_reader = vcf.Reader(f_primary)

    confirming_reader = []
    for f in f_confirming:
        confirming_reader.append(vcf.Reader(f))

    primary_reader.infos["VariantCalledBy"] = _Info(
        "VariantCalledBy", ".", "String", "variant callers that called the variant", "caller ", "0.1"
    )

    # add GT if primary is ST (strelka)
    if primary_caller_name == "ST":
        primary_reader.formats['GT'] = _Format('GT', 1, 'String', 'Genotype')
    if primary_caller_name in ["ST", "VS"]:
        primary_reader.formats['AF'] = _Format('AF', 1, 'Float', 'Allele Frequency')


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

                if primary_caller_name == "VS":  # Varscan case
                    if 'FREQ' in primary_rec.FORMAT:  # Check if FREQ exists
                        formats = primary_rec.FORMAT.split(":")  # Split FORMAT string into a list
                        if "FREQ" in formats:
                            freq_pos= formats.index("FREQ")
                            formats.remove("FREQ")
                        if "AF" not in formats:
                            formats.insert(freq_pos, "AF")  # insert AF after GT

                        primary_rec.FORMAT = ":".join(formats)  # Join the list back into a string

                        for call in primary_rec.samples:
                            if 'FREQ' in call.data._fields:
                                freq_string = call['FREQ']

                                try:
                                    freq = round(float(freq_string[:-1]) / 100, 3)  # Convert percentage string to float
                                    # Update call data
                                    new_fields = list(call.data._fields)
                                    if 'AF' not in call.data._fields:
                                        new_fields.insert(freq_pos, 'AF') # Ensure the field exists before creating the tuple
                                    if 'FREQ' in new_fields:
                                        new_fields.remove('FREQ')

                                    NewCallData = vcf.model.make_calldata_tuple(new_fields)

                                    values = [getattr(call.data, field) for field in call.data._fields if field != 'FREQ']
                                    values.insert(new_fields.index('AF'), freq)
                                    new_call_data = NewCallData(*values)

                                    call.data = new_call_data
                                except (TypeError, ValueError):  # Handle potential errors in conversion
                                    print(f"WARNING: Could not convert FREQ value '{freq_string}' "
                                          f"for sample {call.sample} at {primary_rec.CHROM}:{primary_rec.POS}")
                                    pass  # or handle differently if you don't want to ignore


                # add genotype if primary_caller_name is ST (strelka)
                if primary_caller_name == "ST":
                    if 'GT' not in primary_rec.FORMAT or 'AF' not in primary_rec.FORMAT:
                        parts = primary_rec.FORMAT.split(':')
                        if 'GT' in parts:
                            gt_index = parts.index('GT')
                            parts.insert(gt_index + 1, 'AF')
                            primary_rec.FORMAT = ':'.join(parts)
                        else:  # Handle the case where GT might not exist.
                            primary_rec.FORMAT = 'GT:AF:' + primary_rec.FORMAT

                    # get REF and ALT count keys
                    if not primary_rec.is_indel:
                        ref_base = primary_rec.REF
                        alt_base = primary_rec.ALT[0]  # Assuming only one ALT allele for simplicity
                        ref_counts_key = f"{ref_base}U"
                        alt_counts_key = f"{alt_base}U"
                    else:
                        ref_counts_key = "TAR"
                        alt_counts_key = "TIR"

                    for call in primary_rec.samples:
                        try:
                            ref_counts = call[ref_counts_key]
                            alt_counts = call[alt_counts_key]

                            # Extract tier 1 counts (first comma-delimited value or first element if already a list)
                            tier1_ref_counts = int(ref_counts.split(',')[0]) if isinstance(ref_counts, str) else int(ref_counts[0])
                            tier1_alt_counts = int(alt_counts.split(',')[0]) if isinstance(alt_counts, str) else int(alt_counts[0])

                            # Calculate allele frequency
                            if tier1_ref_counts + tier1_alt_counts > 0:
                                af = tier1_alt_counts / (tier1_ref_counts + tier1_alt_counts)
                                af = round(af, 3)
                            else:
                                af = 0.0

                            if hasattr(call.data, 'AF'):
                                call.data.AF = af
                            else:
                                new_fields = ['AF'] + list(call.data._fields)
                                NewCallData = vcf.model.make_calldata_tuple(new_fields)
                                new_call_data = NewCallData(af, *call.data)
                                call.data = new_call_data

                        except KeyError:
                            # Handle cases where the required FORMAT fields are missing
                            print(f"WARNING: Missing {ref_counts_key} or {alt_counts_key} "
                                f"for sample {call.sample} in record {primary_rec.CHROM}:{primary_rec.POS}")
                            call.data.AF = None

                        # genotype = "0/0" if call.sample == normal_sample_name else "0/1"
                        if call.sample == normal_sample_name:
                            genotype =  primary_rec.INFO.get('NT')
                            if genotype == 'ref':
                                genotype = '0/0'
                            elif genotype == 'het':
                                genotype = '1/0'
                            elif genotype == 'hom':
                                genotype = '1/1'
                            else:
                                genotype = './.'
                        else:
                            genotype = primary_rec.INFO.get('SGT')
                            if genotype:
                                genotype = '0/1'
                            else:
                                genotype = './.'
                        if hasattr(call.data, 'GT'):
                            call.data.GT = genotype
                        else:
                            new_fields = ['GT'] + list(call.data._fields)
                            NewCallData = vcf.model.make_calldata_tuple(new_fields)
                            new_call_data = NewCallData(genotype, *call.data)
                            call.data = new_call_data

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
    parser.add_argument(
        "--normal_sample_name", required=True, type=str, help="Name of the normal sample",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    var_count = make_hc_somatic_vars(
        args.primary, args.primary_name, args.confirming, args.confirming_names, args.out_vcf, args.out_single_vcf
    )

    print("Total confirmed:\t" + str(var_count))
