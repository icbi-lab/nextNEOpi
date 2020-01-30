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
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import argparse
import os
import sys

import vcf
from vcf.parser import _Info
from vcf.utils import walk_together


def make_hc_somatic_vars(fm2, fm1, fvs, fout, fout_single, priority):
    """Picks the given VARs from the priority vcf to a new
    VCF if they are confirmed by one of the others."""

    m2_reader = vcf.Reader(fm2)
    m1_reader = vcf.Reader(fm1)
    vs_reader = vcf.Reader(fvs)

    if priority == 'm2':
        priority_reader = m2_reader
        priority_caller_name = 'M2'
        confirming1_reader = m1_reader
        confirming1_call_name = 'M1'
        confirming2_reader = vs_reader
        confirming2_call_name = 'VS'
    elif priority == 'm1':
        priority_reader = m1_reader
        priority_caller_name = 'M1'
        confirming1_reader = m2_reader
        confirming1_call_name = 'M2'
        confirming2_reader = vs_reader
        confirming2_call_name = 'VS'
    elif priority == 'vs':
        priority_reader = vs_reader
        priority_caller_name = 'VS'
        confirming1_reader = m2_reader
        confirming1_call_name = 'M2'
        confirming2_reader = m1_reader
        confirming2_call_name = 'M1'
    else:
        print("Unknown priority vcf: " + priority)
        sys.exit(1)

    priority_reader.infos['VariantCalledBy'] = _Info('VariantCalledBy', '.', 'String', 'variant callers that called the variant', 'caller ', '0.1')

    # some sanity checks
    if priority_reader.samples != confirming1_reader.samples != confirming2_reader.sample:
        raise ValueError("Input VCF files must have the same sample column "
                "headers.")
    if sorted(priority_reader.contigs.keys()) != sorted(confirming1_reader.contigs.keys()) != sorted(confirming2_reader.contigs.keys()):
        raise ValueError("Input VCF files must denote the same contigs.")

    out_writer = vcf.Writer(fout, priority_reader)
    out_writer_single = vcf.Writer(fout_single, priority_reader)

    for priority_rec, confrirming1_rec, confrirming2_rec in walk_together(priority_reader, confirming1_reader, confirming2_reader):
        confirmed = False
        if priority_rec is not None:
            if confrirming1_rec is not None and confrirming2_rec is not None:
                priority_rec.add_info('VariantCalledBy', priority_caller_name + ',' + confirming1_call_name + ',' + confirming2_call_name)
                confirmed = True
            elif confrirming1_rec is not None and confrirming2_rec is None:
                priority_rec.add_info('VariantCalledBy', priority_caller_name + ',' + confirming1_call_name)
                confirmed = True
            elif confrirming1_rec is None and confrirming2_rec is not None:
                priority_rec.add_info('VariantCalledBy', priority_caller_name + ',' + confirming2_call_name)
                confirmed = True
            else:
                priority_rec.add_info('VariantCalledBy', priority_caller_name)
                out_writer_single.write_record(priority_rec)

            if confirmed:
                out_writer.write_record(priority_rec)

if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(description='Create vcf that contains only vars confirmed by at least one other caller')

    def _file_read(fname, mode='rt'):
        """Returns an open file handle if the given filename exists."""
        if not os.path.exists(fname):
            parser.error("File '{0}' not found.".format(fname))
        _, file_extension = os.path.splitext(fname)
        if file_extension == ".gz":
            mode = 'rb'
        return open(os.path.join(os.path.dirname(__file__), fname), mode)

    def _file_write(fname):
        """Returns an open file handle if the given filename exists."""
        return open(fname, 'w')


    parser.add_argument('--priority', required=True, type=str, choices=['m2', 'm1', 'vs'],
            help='Which caller\'s variant is kept when a variant is called in '
            'both files')
    parser.add_argument('--m2_vcf', required=True, type=_file_read, help='VCF file produced by '
            'Mutect2')
    parser.add_argument('--m1_vcf', required=True, type=_file_read, help='VCF file produced by '
            'Mutect1')
    parser.add_argument('--vs_vcf', required=True, type=_file_read, help='VCF file produced by '
            'Varscan2')
    parser.add_argument('--out_vcf', required=True, type=_file_write, help='VCF file for confirmed vars produced by '
            'this tool')
    parser.add_argument('--out_single_vcf', required=True, type=_file_write, help='VCF file for unconfirmed vars produced by '
            'this tool')

    parser.add_argument('--version', action='version', version='%(prog)s ' +
            __version__)

    args = parser.parse_args()

    make_hc_somatic_vars(args.m2_vcf, args.m1_vcf, args.vs_vcf, args.out_vcf, args.out_single_vcf, args.priority)
