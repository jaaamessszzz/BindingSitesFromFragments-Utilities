#!/usr/bin/env python3

"""
Determine fragment representation in PDB:
    * Number of ligands represented
    * Number of unique ligand-protein complexes
    * Number of unique mappings

Usage:
    fragment_representation <pdb_dir>

Arguments:
    <pdb_dir>
        Directory containing transformed pdb-fragment alignments
"""
import os
import sys
import re
import docopt

def pdb_check(dir, base_only=False):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if path.endswith('.pdb'):
            if base_only:
                yield file
            else:
                yield path

def consolidate_pdb_metrics(pdb_dir):
    """
    Consolidate PDB metrics!
    :return:
    """
    ligand_set = set()
    scaffold_set = set()
    scaffold_ligand_complex_set = set()
    mappings = set()

    for pdb in pdb_check(pdb_dir, base_only=True):

        split_cluster = re.split('-|_', pdb)
        ligand = split_cluster[1]
        scaffold = split_cluster[0]

        ligand_set.add(ligand)
        scaffold_set.add(scaffold)
        scaffold_ligand_complex_set.add(f'{scaffold}_{ligand}')
        mappings.add(pdb)

    print(f'{len(ligand_set)} unique ligands')
    print(f'{len(scaffold_set)} unique scaffolds')
    print(f'{len(scaffold_ligand_complex_set)} unique ligand-scaffold complexes')
    print(f'{len(mappings)} unique mappings onto fragment')

if __name__ == '__main__':

    args = docopt.docopt(__doc__)
    consolidate_pdb_metrics(args['<pdb_dir>'])