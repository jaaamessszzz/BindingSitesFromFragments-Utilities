#!/usr/bin/env python3

"""
Generate part plasmid input sequences from design PDBs in csv format.

Usage:
    sequences_from_pdb dimer <design_PDB_dir> <spec_json> [--output_name]
    sequences_from_pdb monomer <design_PDB_dir> [<spec_json>] [--output_name]

Arguments:
    monomer
        Generate sequences for monomer designs. (Default: generates part 3 inserts)

    dimer
        Generate sequences for dimer designs

    <spec_json>
        JSON defining:
            * part type ( 3a | 3b ) for each chain in the design
            * Reversions/point mutations
        Required for dimer designs, optional for monomers.

    <design_PDB_dir>
        Directory containing design PDBs

Options:
    --output_name -o
        Name for output .csv
"""
import sys
import os
import re

import copy
import json
from pprint import pprint

import docopt
import prody
import pandas as pd

file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(file_dir, '..'))
from utils import pdb_check

# BsaI SITES AND OVERHANGS ONLY!!!
part_prefix_suffix = {'prefix': {'3a': 'GGTCTCATATG',
                                 '3b': 'GGTCTCATTCT'
                                 },
                      'suffix': {'3a': 'GGTTCTTGAGACC',  # Extra GG included for GS linker
                                 '3b': 'GGATCCTGAGACC'   # Extra GG included for GS linker
                                 }
                      }

codon_table_df = pd.read_csv(os.path.join(file_dir, 'YeastCodonTable.csv'))

def make_reversions(hv, design_chain, spec_json):
    """
    Make reversions specified in the specification JSON
    :return:
    """
    chain_start_resnum = hv[design_chain][0].getResnum()
    chain_end_resnum = hv[design_chain][-1].getResnum()

    print(f'Chain {design_chain}', chain_start_resnum, chain_end_resnum)

    chain_sequence_list = [res for res in hv[design_chain].getSequence()]

    if spec_json:
        starting_sequence = copy.deepcopy(chain_sequence_list)

        for reversion in spec_json["reversions"]:
            reversion_resnum = reversion[1]
            if chain_start_resnum <= reversion_resnum <= chain_end_resnum:
                reversion_index_in_sequence = reversion_resnum - chain_start_resnum
                if reversion[0] != 'X':
                    if chain_sequence_list[reversion_index_in_sequence] in reversion[0]:
                        chain_sequence_list[reversion_index_in_sequence] = reversion[2]
                else:
                    chain_sequence_list[reversion_index_in_sequence] = reversion[2]

        # VERIFICATION
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment

        alignments = pairwise2.align.globalms(''.join(starting_sequence), ''.join(chain_sequence_list), 2, -1, -10,
                                              -.5)  # No Gaps!
        print(format_alignment(*alignments[0]))

    return chain_sequence_list

def generate_sequences_from_pdb(design_dir, spec_json=None):
    """
    Generate DNA sequence from input PDB
    :param design_dir:
    :param spec_json:
    :return:
    """
    # --- Load specifications if present --- #
    if spec_json:
        print('Loading specifications...')
        spec_json = json.load(open(args['<spec_json>'], 'r'))

        print(f'The following specifications were found:')
        pprint(spec_json)

    design_dict_list = []

    for index, design in enumerate(pdb_check(design_dir)):

        hv = prody.HierView(prody.parsePDB(design))
        design_chains = set([chain.getChid() for chain in hv.iterChains()]) - set('X')

        current_design = re.split('/|\.', design)[1]
        print(current_design)

        for design_chain in design_chains:
            chain_sequence_list = make_reversions(hv, design_chain, spec_json)

            design_dict = {'design': current_design,
                           'index': index + 1,
                           'chain': design_chain,
                           'part': spec_json['part_types'][design_chain],
                           'AA_sequence': ''.join(chain_sequence_list)
                           }

            design_dict_list.append(design_dict)

    # --- Dump sequences into csv --- #
    df = pd.DataFrame(design_dict_list)
    df.to_csv('20181001-4YDY_Designs-AA_Sequences.csv', index=False)

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    generate_sequences_from_pdb(args['<design_PDB_dir>'], spec_json=args['<spec_json>'])
