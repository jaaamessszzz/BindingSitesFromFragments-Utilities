#!/usr/bin/env python3

"""
Generates a constraint file based on an input binding motif and constraint json.

Usage:
    Prep-cst_from_motif <motif_pdb> <cst_json>

Arguments:
    <motif_pdb>
        PDB of binding motif
        
    <cst_json>
        JSON containing all constraint blocks used for matching
"""

import os
import sys
import subprocess
import json
import re

import docopt
import prody

if __name__ == '__main__':
    
    args = docopt.docopt(__doc__)
    motif_bpdb_name = os.path.basename(args['<motif_pdb>'])
    cst_name = '{0}.cst'.format(motif_bpdb_name.split('.')[0])
    conformer = motif_bpdb_name.split('-')[0]

    with open(args['<cst_json>'], 'r') as cst_json:

        # Get required blocks from motif pdb file name
        cst_json_json = json.load(cst_json)
        relevant_constraint_blocks = re.split('-|\.', motif_bpdb_name)[1].split('_')[1:]

        with open(cst_name, 'w') as cst_file:
            for cst_block_index in relevant_constraint_blocks:
                cst_file.write('CST::BEGIN\n')
                cst_file.write(cst_json_json[conformer][cst_block_index])
                cst_file.write('\n')
                cst_file.write('CST::END\n')

        print('\nGenerated {0}!\n'.format(cst_name))