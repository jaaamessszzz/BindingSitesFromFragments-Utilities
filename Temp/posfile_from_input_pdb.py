#! /usr/bin/env python3

"""
Spits out a posfile for a given pdb. Literally prints all resnums from an input pdb to a text file, using this for biasing
Gurobi solutions for preselected residues.

Usage:
posfile_from_input_pdb <input_pdb>

Arguments:
    <input_pdb>
        Path to input pdb
"""

import os
import docopt
import prody

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    prody_pdb = prody.parsePDB(args['<input_pdb>']).getHierView()
    posfile_name = os.path.basename(os.path.normpath(args['<input_pdb>'])).split('.')[0] + '.pos'

    with open(posfile_name, 'w') as posfile:
        for residue in prody_pdb.iterResidues():
            posfile.write('{}\n'.format(residue.getResnum()))