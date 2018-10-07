'''
PcrGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
import os
import shutil
import sys

from Bio import Seq
from synbiochem.utils import ice_utils, seq_utils

import pandas as pd


def get_ice_files(url, username, password,
                  ice_ids_filename,
                  for_primer, rev_primer,
                  dir_name):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)
    ice_ids = _get_ice_ids(ice_ids_filename)

    seqs_offsets = [pcr(ice_client.get_ice_entry(ice_id).get_seq(),
                        for_primer, rev_primer)
                    for ice_id in ice_ids]

    seqs, _ = zip(*seqs_offsets)

    _mkdirs(dir_name)

    for ice_id, seq in zip(ice_ids, seqs):
        seq_utils.write_fasta({ice_id: seq}, os.path.join(dir_name,
                                                          ice_id + '.fasta'))

    # Get Genbank files for subsequent data analysis:
    for ice_id in ice_ids:
        gb_filename = os.path.join(dir_name, ice_id + '.gb')
        ice_client.get_genbank(ice_id, gb_filename)

    return ice_ids, [len(seq) for seq in seqs]


def pcr(seq, forward_primer, reverse_primer):
    '''Apply in silico PCR.'''
    for_primer_pos = seq.find(forward_primer.upper())

    rev_primer_pos = \
        seq.find(str(Seq.Seq(reverse_primer).reverse_complement().upper()))

    if for_primer_pos > -1 and rev_primer_pos > -1:
        seq = seq[for_primer_pos:] + \
            seq[:rev_primer_pos + len(reverse_primer)]
    elif for_primer_pos > -1:
        seq = seq[for_primer_pos:]
    elif rev_primer_pos > -1:
        seq = seq[:rev_primer_pos + len(reverse_primer)]

    return seq, for_primer_pos


def _get_ice_ids(ice_ids_filename):
    '''Get ICE ids.'''
    with open(ice_ids_filename, 'rU') as ice_ids_file:
        return [line.strip() for line in ice_ids_file]


def _mkdirs(dir_name):
    '''Make directories.'''
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)

    os.makedirs(dir_name)


def main(args):
    '''main method.'''
    ice_ids, lngths = get_ice_files(*args[:-1])

    df = pd.DataFrame()
    df['ice_id'] = ice_ids
    df['length'] = lngths

    df.to_csv(os.path.join(args[-2], args[-1]), index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
