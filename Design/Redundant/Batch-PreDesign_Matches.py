#!/usr/bin/python
"""
Script for submitting jobs to predesign matches for CoupledMoves and BackrubEnsemble
"""

from __future__ import print_function
import os
import subprocess
import shutil
import fileinput

cwd = os.getcwd()

def move_file_to_design_dir(file_to_move, design_dir):
    """
    Moves `file_to_move` into `design_dir`

    :param file_to_move: file to move
    :param design_dir: destination directory for `file_to_move`
    """
    src = os.path.join(os.getcwd(), file_to_move)
    dest = os.path.join(design_dir, file_to_move)
    shutil.copy(src, dest)


# --- Positional Arguments --- #

# Do design first with FixBB, filter, then use the best 5-10 then 10 backrub as starting points

# --- Fixed variables --- #
# Get compound ID
cwd_list = os.getcwd().split('/')
compound_id = cwd_list[cwd_list.index('Compounds') + 1]

compound_dir = '/netapp/home/james.lucas/BindingSitesFromFragments/Compounds/{0}'.format(compound_id)
params_dir = os.path.join(compound_dir, 'params')

# --- Submit jobs for each input PDB --- #
for pdb in os.listdir(cwd):
    if pdb.endswith('.pdb'):

        current_xml = 'BackrubEnsemble-Ensemble.xml'

        # input_pdb = sys.argv[1]
        # params_file_path = sys.argv[2]
        # design_json_path = sys.argv[3]
        # designs_per_task = sys.argv[4]

        # Make/check directories
        design_path = os.path.join(cwd, current_xml.split('.')[0], pdb.split('.')[0])

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

        # Move Design_Template.xml into design directory
        move_file_to_design_dir('Design_Template.xml', design_path)

        # Replace "%%design_xml%%" in Design_Template.xml with correct XML
        template_xml = fileinput.FileInput(os.path.join(design_path, 'Design_Template.xml'), inplace=True)
        for line in template_xml:
            print(line.replace('%%design_xml%%', current_xml), end='')
        template_xml.close()

        # Move submission script into design directory
        current_submit_script = 'Submit-Preparation-{0}.py'.format(current_xml.split('.')[0])
        move_file_to_design_dir(current_submit_script, design_path)

        # Get correct params file from PDB name
        params_file = pdb.split('-')[0][-8:]
        print(params_file)

        # todo: add check for relaxed PDB inputs (CoupledMoves and BackrubEnsemble)? Or always generate ab initio?
        # todo: create dict for appropriate args for each design XML

        args = ['qsub',
                current_submit_script,
                os.path.join(compound_dir, 'Design', pdb),
                os.path.join(params_dir, '{0}.params'.format(params_file)),
                '../../{0}-design_inputs.json'.format(pdb.split('.')[0]),
                '5' # Designs per task
                ]
        print(args)

        process = subprocess.Popen(args, cwd=design_path)
        process.wait()

