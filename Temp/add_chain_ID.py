#!/usr/bin/env python3

"""
Adds chain ID to input PDBs... necessary for adding constraints to Baker Lab monomer scaffolds...

Usage:
    add_chain_ID <Input_PDB> <Chain_ID>

Arguments:
    <Input_PDB>
        Path to input PDB

    <Chain_ID>
        Letter to use for chain ID (CANNOT BE X)
"""
import os
import fileinput

import docopt

# Adapted from https://stackoverflow.com/a/41752999
def replace_str_index(text,index=0,replacement=''):
    return '{0}{1}{2}'.format(text[:index],replacement,text[index+1:])


args = docopt.docopt(__doc__)

pdb_path = args['<Input_PDB>']
chid = args['<Chain_ID>']

# Replace "%%design_xml%%" in Design_Template.xml with correct XML
shit_pdb = fileinput.FileInput(pdb_path, inplace=True)
for line in shit_pdb:
    if line.startswith('ATOM'):
        print(replace_str_index(line, index=21, replacement=chid), end='')
    else:
        print(line, end='')

shit_pdb.close()