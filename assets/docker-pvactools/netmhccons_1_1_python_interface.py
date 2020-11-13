import os
from subprocess import Popen, PIPE
import tempfile

from allele_info.allele_info import is_user_defined_allele, MHCIAlleleData
from iedbtools_utilities.sequence_io import SequenceOutput
from iedbtools_utilities.sequence_manipulation import split_sequence
from logging import getLogger
logger = getLogger(__name__)

EXECUTABLE_NAME = 'netMHCcons'
EXECUTABLE_DIR_PATH = os.path.normpath(os.path.join(os.path.abspath(__file__), os.path.pardir))
EXECUTABLE_FULL_PATH = os.path.join(EXECUTABLE_DIR_PATH, EXECUTABLE_NAME)

def predict_sequence(sequence, allele_length_pair):
    '''Given one protein sequence, break it up into peptides, return their predicted binding scores.'''
    allele, binding_length = allele_length_pair
    peptide_list = split_sequence(sequence, binding_length)
    scores = predict_peptide_list(peptide_list, allele_length_pair)
    return scores

def parse_pickpocket_output(content):
    '''Given pickpocket output, returns a list of scores'''
    scores = []
    for lines in content.split('\n'):
        if 'PEPLIST' in lines:
            data_list = lines.split()
            if data_list[0].isdigit():
                scores.append(float(data_list[5]))
    return scores

def predict_peptide_list(peptide_list, allele_length_pair):
    '''This routine can be directly called so that you do not make a file for each prediction.'''
    infile = tempfile.NamedTemporaryFile(prefix='netmhccons_', suffix='_input', delete=False)
    for peptide in peptide_list:
        infile.write(peptide + "\n")
    infile.close()
    allele, binding_length = allele_length_pair
    
    # Temporary fix.
    stripped_allele_name = strip_allele_name(allele)

    user_defined_allele = is_user_defined_allele(allele)

    if user_defined_allele:
        fasta_allele = SequenceOutput(allele).to_fasta()
        usermhcfile = tempfile.NamedTemporaryFile(
                          prefix='netmhccons_', suffix='_usermhc', delete=False)
        usermhcfile.write(fasta_allele)
        usermhcfile.close()

    if user_defined_allele:
        cmd = [
            EXECUTABLE_FULL_PATH, '-hlaseq', usermhcfile.name, '-length', str(binding_length),
            '-inptype', '1', '-f', infile.name
        ]
    else:
        cmd = [
            EXECUTABLE_FULL_PATH, '-a', stripped_allele_name, '-length', str(binding_length),
            '-inptype', '1', '-f', infile.name
        ]

    logger.info('Calling netmhccons executable:\n{}'.format(' '.join(cmd)))
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    #process_status_code = p.wait()
    analysis_results, ignored_stderr = p.communicate()
    process_status_code = p.returncode
    #scores = parse_pickpocket_output(analysis_results)
#     process_status_code = f.close()
    if process_status_code != 0:
        msg = 'Error calling netmhccons executable:\n{}'.format(' '.join(cmd))
        logger.error(msg)
        raise Exception(msg)
    scores = parse_pickpocket_output(analysis_results)

    os.remove(infile.name)
    if user_defined_allele:
        os.remove(usermhcfile.name)

    if (len(peptide_list) != len(scores)):
        msg = "len(peptide_list) != len(scores) -- {} != {}".format(len(peptide_list), len(scores))
        logger.error(msg + '\n' + analysis_results)
        raise Exception(msg)

    if cmd == None:  # According to python doc, return status of '0' is actually 'None' using os.popen.
        msg = "netMHCconsPredictor did not execute. Command used:\n{}".format(cmd)
        raise Exception(msg)

    return tuple(scores)

def strip_allele_name(allele_name):
    """ | *brief*: Temporary hack to get the allele name right for netmhccons executable.
        | *author*: Dorjee
        | *created*: 2016-09-13

        TODO: A more permanent solution would be to create a column in the database for canonical allele name.
    """
    miad = MHCIAlleleData()
    species = miad.get_species_for_allele_name(allele_name=allele_name)
    if species in ['macaque', 'pig']:
        stripped_allele_name = allele_name.replace('*',':')
    else:
        stripped_allele_name = allele_name.replace('*', '')
    return stripped_allele_name
