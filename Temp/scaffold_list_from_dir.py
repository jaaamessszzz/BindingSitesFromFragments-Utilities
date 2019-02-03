#!/usr/bin/env python3

"""
Make a text file that lists all PDBs to use for matching

Usage:
    scaffold_list_from_dir <dir_name> <scaffold_file_name>

Arguments:

    <dir_name>
        Directory containing scaffolds

    <scaffold_file_name>
        File to save scaffold list
"""
import os
import sys

import docopt

def main():
    args = docopt.docopt(__doc__)

    with open(args['<scaffold_file_name>'], 'w') as file:
        for pdb in sorted(os.listdir(args['<dir_name>'])):
            if '.pdb' in pdb:
                file.write('{0}\n'.format(pdb))

if __name__ == '__main__':
    main()