#! /usr/bin/python3

"""
Generate matcher constraints for an existing protein-ligand complex

Usage:
  constraints walk <root_dir> [options]
  constraints <pdb_path> <ligand_code> <ligand_position> <positions>... [options]

Arguments:
  <pdb_path>            Path to PDB for scoring
  <positions>           Positions of residues in input PDB to generate constraints
  <ligand_code>         Chemical component identifier for ligand
  <ligand_position>     Position of ligand to calculate interface energies

Options:
  --clean, -c                                   Clean input PDB before assessing energies
  -t=<tolerance>, --tolerance=<tolerance>       Matcher tolerance (in degrees) per constraint
  -s=<samples>, --samples=<samples>             Matcher samples per constraint
  -g, --greasy                                  Greasy sampling
  -p=<prefix>, --params=<prefix>                Use <prefix> for writing constraints in place of ligand_code

"""

import docopt
import io
import os
import shutil
import pandas as pd
import prody

from BindingSitesFromFragments.motifs import Generate_Constraints

if __name__ == '__main__':

    args = docopt.docopt(__doc__)

    tolerance = int(args['--tolerance']) if args['--tolerance'] else 5
    samples = int(args['--samples']) if args['--samples'] else 1
    greasy = args['--greasy']

    params_dir = 'params'
    os.makedirs(params_dir, exist_ok=True)

    if not args['walk']:
        pdb_path = args['<pdb_path>']
        positions_list = [int(position) for position in args['<positions>']]
        ligand_code = args['<ligand_code>']
        ligand_position = int(args['<ligand_position>'])


        assert ligand_position not in positions_list

        constraint_generator = Generate_Constraints(ligand_code)
        pdb_prody = prody.parsePDB(pdb_path)
        ligand_prody = pdb_prody.select(f'resname {ligand_code} and resnum {ligand_position}').copy()
        constraint_prefix = args['--params'] if args['--params'] else ligand_code

        with open(f'{constraint_prefix}-existing-{"_".join(args["<positions>"])}.cst', 'w') as current_constraint_file:

            for position in positions_list:
                residue_prody = pdb_prody.select(f'resnum {position}').copy()
                constraint_atoms_dict, constraint_block = constraint_generator.generate_single_constraint_block_base(residue_prody, ligand_prody,
                                                                                                                     angle_A_tolerance=tolerance,
                                                                                                                     angle_B_tolerance=tolerance,
                                                                                                                     torsion_A_tolerance=tolerance,
                                                                                                                     torsion_AB_tolerance=tolerance,
                                                                                                                     torsion_B_tolerance=tolerance,
                                                                                                                     torsion_constraint_sample_number=samples,
                                                                                                                     angle_constraint_sample_number=samples,
                                                                                                                     greasy_sampling=greasy)

                current_constraint_file.write('CST::BEGIN\n')
                current_constraint_file.write('\n'.join(constraint_block))
                current_constraint_file.write('\n')
                current_constraint_file.write('CST::END\n')

    else:
        for root, dirs, files in os.walk(args['<root_dir>']):
            print(root)
            params_files = [file for file in files if file.endswith('.params')]
            for params in params_files:
                shutil.copy(os.path.join(root, params), params_dir)
                param_name, _ = os.path.splitext(params)
                param_split = param_name.split('_')
                lig_id = None

                if len(param_split) > 2:
                    ligand = param_split[0]
                    complex = param_split[1]
                    lig_id = param_split[2]
                    interaction_df = pd.read_csv(os.path.join(root, f'{complex}_{lig_id}-LigandInterfaceEnergy.csv'))
                    pdb_prody_path = os.path.join(root, f'{complex}_{lig_id}-clean.pdb')

                else:
                    ligand = param_split[0]
                    complex = param_split[1]
                    interaction_df = pd.read_csv(os.path.join(root, f'{complex}-LigandInterfaceEnergy.csv'))
                    pdb_prody_path = os.path.join(root, f'{complex}-clean.pdb')

                positions_list = interaction_df.head(3)['position']
                print(interaction_df.head(3))
                print(ligand, complex, lig_id)

                constraint_generator = Generate_Constraints(ligand)
                pdb_prody = prody.parsePDB(pdb_prody_path)
                ligand_position = max(pdb_prody.getResnums())

                print(ligand, ligand_position)
                ligand_prody = pdb_prody.select(f'resname {ligand} and resnum {ligand_position}').copy()

                with open(f'{param_name}-existing-{"_".join([str(a) for a in positions_list])}.cst', 'w') as current_constraint_file:

                    for position in positions_list:
                        residue_prody = pdb_prody.select(f'resnum {position}').copy()
                        constraint_atoms_dict, constraint_block = constraint_generator.generate_single_constraint_block_base(
                            residue_prody, ligand_prody,
                            angle_A_tolerance=tolerance,
                            angle_B_tolerance=tolerance,
                            torsion_A_tolerance=tolerance,
                            torsion_AB_tolerance=tolerance,
                            torsion_B_tolerance=tolerance,
                            torsion_constraint_sample_number=samples,
                            angle_constraint_sample_number=samples,
                            greasy_sampling=greasy)

                        current_constraint_file.write('CST::BEGIN\n')
                        current_constraint_file.write('\n'.join(constraint_block))
                        current_constraint_file.write('\n')
                        current_constraint_file.write('CST::END\n')
