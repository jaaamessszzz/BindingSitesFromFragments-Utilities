#!/usr/bin/python
"""
Script for submitting design strategy simulations for all motif-scaffold matches
Run when ssh to iqint

Start with directory containing this script and all scaffolds and associated design JSON
Each Submit script requires positional arguments for match PDB, params file path, and design JSON

+-- Compounds
    +-- REN
        +-- Design
            +-- Submit-DesignFilterCutoffs.py
            +-- <All Design XMLs>
            +-- <All Submission Scripts>

            # This script will use names of PDBs in this directory to find inputs for all design strategies
            +-- <Raw_input_dir>
                +-- Match_1.pdb
                +-- Match_2.pdb

            # PDBs in this directory will be used for BackrubEnsemble design
            +-- Inputs-BackrubEnsemble
                +-- PDB_1.pdb
                +-- PDB_2.pdb

            # PDBs in this directory will be used for CoupledMoves design
            +-- Inputs-Predesigned
                +-- PDB_1.pdb
                +-- PDB_2.pdb

            # These directories are created by Submit-DesignFilterCutoffs.py
            # Jobs submitted from these directories
            +-- <Matched PDB Name>
                +-- <Design Strategy>

                    # XML is edited to replace "%%design_xml%%" with the desired design XML
                    +-- Design_Template.xml

        +-- params
            +-- REN_00XX.params

From the Design directory, run: scl enable python27 'Batch-DesignFilterCutoffs.py <Raw_input_dir>

This script will use the names of the PDBs in <Raw_input_dir> to determine which inputs should be used for all design
strategies

"""
from __future__ import print_function
import os
import sys
import re
import subprocess
import shutil
import fileinput

def move_file_to_design_dir(file_to_move, design_dir):
    """
    Moves `file_to_move` into `design_dir`

    :param file_to_move: file to move
    :param design_dir: destination directory for `file_to_move`
    """
    src = os.path.join(os.getcwd(), file_to_move)
    dest = os.path.join(design_dir, file_to_move)
    shutil.copy(src, dest)

# todo: CoupledMoves and BackrunEnsemble require FastRelaxed inputs...
design_strategies = ['FastDesign_HBNet.xml',
                     'FastDesign_Only.xml',
                     # 'CoupledMoves.xml',
                     # 'BackrubEnsemble.xml',
                     ]

# --- Positional Arguments --- #
design_inputs = sys.argv[1]
foldtree_path = sys.argv[2]
cstfile_path = sys.argv[3]

# --- Fixed variables --- #
# Get compound ID
cwd = os.getcwd()
cwd_list = cwd.split('/')
compound_id = cwd_list[cwd_list.index('Compounds') + 1]

# Directory with params files
params_dir = '/netapp/home/james.lucas/BindingSitesFromFragments/Compounds/{0}/params'.format(compound_id)

# --- Submit jobs for each input PDB --- #
for pdb in os.listdir(design_inputs):
    if pdb.endswith('.pdb'):

        current_pdb_name = pdb.split('.')[0]

        # Create Design Directory
        try:
            os.makedirs(current_pdb_name)
        except Exception as e:
            print(e)
            print("WARNING: Directory for {0} already exists! Continuing anyways...\n\n".format(current_pdb_name))

        # Submit jobs for each design strategy
        for current_xml in design_strategies:

            # Make/check directories
            design_path = os.path.join(cwd, current_pdb_name, current_xml.split('.')[0])

            # FastDesign Methods can have everything made from scratch using matcher output
            if current_xml in ['FastDesign_HBNet.xml', 'FastDesign_Only.xml']:

                # Only continue if directory doesn't exist, else continue
                # The logic being that if the directory exists, something was already submitted
                # Else, just delete/rename the existing directory and run this again
                if not os.path.exists(design_path):
                    os.makedirs(design_path)
                else:
                    print('\n\n'
                          'Skipping {0} - directory exists!\n'
                          'Either delete or rename the existing directory and resubmit!'
                          '\n\n'.format(current_xml))
                    continue

                qsub_pdb_argument = os.path.join(cwd, design_inputs, pdb)

            # CoupledMoves and BackrubEnsemble require pre-generated FixBB outputs to alleviate clashes introduced by
            # matched motif residues
            elif current_xml == 'CoupledMoves.xml':

                # If path with pre-generated inputs does not exist, move on
                if not os.path.exists(os.path.join(cwd, 'Inputs-Predesigned', current_pdb_name)):
                    print('\n\n'
                          'Skipping {0} - required input directory does not exist!'.format(current_xml))
                    continue

                # todo: get rid of repeat code...
                if not os.path.exists(design_path):
                    os.makedirs(design_path)
                else:
                    print('\n\n'
                          'Skipping {0} - directory exists!\n'
                          'Either delete or rename the existing directory and resubmit!'
                          '\n\n'.format(current_xml))
                    continue

                qsub_pdb_argument = os.path.join(cwd, 'Inputs-Predesigned', current_pdb_name)
                designs_per_task = '5'

            elif current_xml == 'BackrubEnsemble.xml':
                # If path with pre-generated inputs does not exist, move on
                if not os.path.exists(os.path.join(cwd, 'Inputs-BackrubEnsemble', current_pdb_name)):
                    print('\n\n'
                          'Skipping {0} - required input directory does not exist!'.format(current_xml))
                    continue

                # todo: get rid of repeat code...
                if not os.path.exists(design_path):
                    os.makedirs(design_path)
                else:
                    print('\n\n'
                          'Skipping {0} - directory exists!\n'
                          'Either delete or rename the existing directory and resubmit!'
                          '\n\n'.format(current_xml))
                    continue

                qsub_pdb_argument = os.path.join(cwd, 'Inputs-BackrubEnsemble', current_pdb_name)
                designs_per_task = '5'

            # --- Move relevant Rosettascript XMLs to design directory --- #

            # Move Design_Template.xml into design directory
            move_file_to_design_dir('Design_Template.xml', design_path)

            # Special case for MC_HBNet since I didn't plan this well...
            if current_xml == 'FastDesign_HBNet.xml':
                move_file_to_design_dir('FastDesign_HBNet.xml', design_path)

            xml_substitution = 'FastDesign_Only.xml' if current_xml == 'FastDesign_HBNet.xml' else current_xml

            # Replace "%%design_xml%%" in Design_Template.xml with correct XML
            template_xml = fileinput.FileInput(os.path.join(design_path, 'Design_Template.xml'), inplace=True)
            for line in template_xml:
                line = line.replace('%%design_xml%%', xml_substitution)
                print(line, end='')

            template_xml.close()

            # Move submission script into design directory
            current_submit_script = 'Submit-Design-{0}.py'.format(current_xml.split('.')[0])
            move_file_to_design_dir(current_submit_script, design_path)

            # Get correct params file from PDB name
            params_file = pdb.split('-')[0][-8:]

            # Get correct constraint file name from PDB name
            possible_cstfiles = [cst for cst in os.listdir(cstfile_path) if cst.endswith('.cst')]
            lmao = [cst.split('.')[0] in pdb for cst in possible_cstfiles]
            cstfile_name = possible_cstfiles[lmao.index(True)]

            # todo: create dict for approproate args for each design XML

            # --- Submit things --- #

            args = ['qsub',
                    current_submit_script,
                    qsub_pdb_argument,
                    os.path.join(params_dir, '{0}.params'.format(params_file)),
                    os.path.join(cwd, '{0}-design_inputs.json'.format(pdb.split('.')[0])),
                    os.path.join(foldtree_path, '{0}-FoldTree.ft'.format(pdb.split('.')[0])),
                    os.path.join(cstfile_path, '{0}'.format(cstfile_name))
                    ]

            if current_xml in ['BackrubEnsemble.xml', 'CoupledMoves.xml']:
                args.append(designs_per_task)

            print('\n')
            print(args)

            process = subprocess.Popen(args, cwd=design_path)
            process.wait()