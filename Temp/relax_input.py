#! /usr/bin/python3

"""
Relax an input PDB File

Usage:
  relax <pdb_path> [options]

Arguments:
  <pdb_path>            Path to PDB to relax

Options:
  -s, --start_coords    Relax to start coords


"""

import docopt
import io
import os
import pandas as pd
import prody

import pyrosetta
from pyrosetta import rosetta

if __name__ == '__main__':

    args = docopt.docopt(__doc__)
    input_pdb_path = args['<pdb_path>']

    sfxn = rosetta.core.scoring.get_score_function()

    relax = pyrosetta.rosetta.protocols.relax.FastRelax(sfxn, 'MonomerRelax2019')
    relax.constrain_relax_to_start_coords(args['--start_coords'])
    relax.apply(input_pdb_path)