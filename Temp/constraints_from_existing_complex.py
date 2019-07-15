#! /usr/bin/python3

"""
Generate matcher constraints for an existing protein-ligand complex

Usage:
  constraints <pdb_path> <positions>... <ligand_code> <ligand_position> [options]

Arguments:
  <pdb_path>            Path to PDB for scoring
  <positions>           Positions of residues in input PDB to generate constraints
  <ligand_code>         Chemical component identifier for ligand
  <ligand_position>     Position of ligand to calculate interface energies

Options:
  --clean, -c           Clean input PDB before assessing energies

"""

import docopt
import io
import os
import pandas as pd
import prody

from BindingSitesFromFragments.motifs import Generate_Constraints

if __name__ == '__main__':

    args = docopt.docopt(__doc__)

    pdb_path = args['<pdb_path>']
    positions_list = [int(position) for position in args['<positions>']]
    ligand_code = args['<ligand_code']
    ligand_position = int(args['<ligand_position>'])

    assert ligand_position not in positions_list

    constraint_generator = Generate_Constraints(ligand_code)
    pdb_prody = prody.parsePDB(pdb_path)
    ligand_prody = pdb_prody.select(f'resname {ligand_code} and resnum {ligand_position}')

    with open(f'{ligand_code}-{"_".join(positions_list)}.cst', 'w') as current_constraint_file:

        for position in positions_list:
            residue_prody = pdb_prody.select(f'resnum {position}').copy()
            constraint_atoms_dict, constraint_block = constraint_generator.generate_single_constraint_block_base(residue_prody, ligand_prody)

            current_constraint_file.write('CST::BEGIN\n')
            current_constraint_file.write('\n'.join(constraint_block))
            current_constraint_file.write('\n')
            current_constraint_file.write('CST::END\n')