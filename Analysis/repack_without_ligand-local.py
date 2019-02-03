#!/usr/bin/python3

import sys
import subprocess
import os
import shutil
import json
import re

import docopt

def determine_matched_residue_positions(match_pdb_path):
    """
    Parse the filename of the match PDB to determine IDs and positions of match residues
    :return:
    """
    positions_block = os.path.basename(os.path.normpath(match_pdb_path)).split('_')[2]
    resnames = [a for a in re.split("[0-9]*", positions_block) if a]
    resnums = [int(a) for a in re.split("[a-zA-Z]*", positions_block) if a]

    return [(a, b) for a, b in zip(resnames, resnums)]

########################
# Positional arguments #
########################

input_pdb = sys.argv[1]
params_file_path = sys.argv[2]
design_json_path = sys.argv[3]
foldtree_path = sys.argv[4]

input_pdb_base = os.path.basename(os.path.normpath(input_pdb))
design_json = json.load(open(design_json_path, 'r'))

##########################
# Start submitting tasks #
##########################


print(input_pdb_base)
print(determine_matched_residue_positions(input_pdb_base))
print(params_file_path)
print(design_json)
print(','.join([str(a) for a in design_json['design_residue_list']]))

arg = ['/netapp/home/james.lucas/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease',
       '-database',
       '/netapp/home/james.lucas/Rosetta/main/database',
       '-s',
       input_pdb,
       '-extra_res_fa',
       params_file_path,
       '-parser:protocol',
       'Repack_without_ligand.xml',
       '-run:preserve_header',  # Required (?) for applying matcher constraints
       # '-use_input_sc',
       '-flip_HNQ',
       '-no_optH',
       'false',
       '-out:prefix',
       '{0}-'.format(sge_task_id + 1),
       '-nstruct',
       '5',
       '-parser:script_vars',
       'motif_residues={0}'.format(','.join([str(resnum[1]) for resnum in determine_matched_residue_positions(input_pdb_base)])),
       'design_positions={0}'.format(','.join([str(a) for a in design_json['design_residue_list']])),
       # 'design_xml={0}'.format('FastDesign.xml')
       'foldtree_file={0}'.format(foldtree_path)
       ]

print(' '.join(arg))

rosetta_process = subprocess.Popen(arg, cwd=os.getcwd())
return_code = rosetta_process.wait()

print('Task return code:', return_code, '\n')

