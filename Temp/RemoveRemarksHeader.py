#!/usr/bin/env python3

"""
Remove prody selection remarks from PDB headers

Usage:
    RemoveRemarksHeader <pdb_dir>

Arguments:
    <pdb_dir>           Directory containing scaffolds to clean
"""

import os
import io
import docopt

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    pdb_dir = args['<pdb_dir>']

    for pdb in os.listdir(pdb_dir):
        if pdb.endswith('.pdb'):

            clean_pdb = io.StringIO()
            with open(os.path.join(pdb_dir, pdb), 'r') as pdb_stream:
                for line in pdb_stream:
                    if line.startswith('REMARK'):
                        continue
                    clean_pdb.write(line)

            with open(os.path.join(pdb_dir, pdb), 'w') as pdb_stream:
                a = clean_pdb.getvalue()
                print(a)
                pdb_stream.write(a)

