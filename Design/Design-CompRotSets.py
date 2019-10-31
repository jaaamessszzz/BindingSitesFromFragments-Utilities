#! /netapp/home/james.lucas/anaconda3/bin/python3
#$ -S /netapp/home/james.lucas/anaconda3/bin/python3
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-60
#$ -j y
#$ -l mem_free=15G
#$ -l scratch=10G

import os
import sys
import tempfile
import shutil

from docopt import docopt
from BindingSitesFromFragments.design import fuzzball_composition_design

def design(args):
    """
    Design the context around the ligand in a protein-ligand complex.
    This is a modified version of the BSFF design command adapted for benchmarking runs on the QB3 Cluster.

    Usage: bsff design <ligand_conformer_path> <match_path> <match_residue_map> <design_config_json> <params_path> [options]

    Arguments:
      <ligand_conformer_path>           Path to ligand PDB generated by molfile_to_params.py
      <match_path>                      Path to scaffold-ligand match PDB
      <match_residue_map>               Pickle file generated by a bsff.assemble iteration
      <design_config_json>              Configuration JSON generated by the Prep-HBNet_Design.py utility
      <params_path>                     Path to the ligand params file generated by molfile_to_params.py

    Options:
      -d=<dir>, --design_dir=<dir>                  Directory to dump designs into
      -l=<limit>, --rotset_limit=<limit>            Limit of rotamers added per residue identity per position
      -m, --apply_minimization                      Apply minimization before accepting/rejecting rotamers
      -n=<nstruct>, --nstruct=<nstruct>             Number of design trajectories to attempt
      -w=<weight>, --special_rot_weight=<weight>    Weight for complementary rotamers during design
      -x, --disable_rotamersets                     Don't use complementary rotamersets during design
    """
    ligand_conformer_path = args['<ligand_conformer_path>']
    match_path = args['<match_path>']
    match_residue_map = args['<match_residue_map>']
    design_config_json = args['<design_config_json>']
    params_path = args['<params_path>']

    # nstruct 50
    # 250 designs/weight/scaffold
    # weights: 12 total, arange(-10,0,2) U arange(-4, -2, 0.25)

    weights = [-10, -8, -6, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0]
    taskid = int(os.environ.get('SGE_TASK_ID'))
    total_tasks = int(os.environ.get('SGE_TASK_LAST'))
    special_rot_weight = weights[(taskid - 1) // int(total_tasks/len(weights))]

    design_dir = args['--design_dir'] if args['--design_dir'] else f'Designs-taskid_{taskid}-weight_{special_rot_weight}'
    rotset_limit = int(args['--rotset_limit']) if args['--rotset_limit'] else 50
    nstruct = int(args['--nstruct']) if args['--nstruct'] else 50

    if special_rot_weight == 0:
        use_complementary_rotsets = False

    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as working_temp_dir:
        design_dir = os.path.join(working_temp_dir, design_dir)
        fuzzball_composition_design(ligand_conformer_path, match_path, match_residue_map, design_config_json, params_path,
                                    rotset_limit=rotset_limit, nstruct=nstruct, apply_minimization=args['--apply_minimization'],
                                    special_rot_weight=special_rot_weight, designdir=design_dir, use_complementary_rotsets=use_complementary_rotsets)

        shutil.move(design_dir, os.environ.get('SGE_O_WORKDIR'))

if __name__ == '__main__':
    argv = sys.argv[1:]
    design(docopt(design.__doc__, argv=argv))