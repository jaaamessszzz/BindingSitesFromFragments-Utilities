#!/usr/bin/env python3
"""
Prepares inputs for the next round of match iteration. This script assumes that a directory containing Pymol sessions of
quality matches named "[match_name]-pretty.pse" exists.

This script will generate the following:

    * New directory with copies of quality match .pdbs from <Match_Dir>
    * .csv and .txt files required for BSFF matching iterations

Idea: generate new posfiles for iterations where only positions directly adjacent to the ligand position are considered.
This should cut down a lot of wasted time considering match positions we know will not be productive.

Usage:
    prepare_match_iteration <Pymol_Dir> <Match_Dir> [--monomer] [--posfile <ligand_ID> <scaffold_file_path>]

Arguments:
    <Pymol_Dir>
        Path to directory containing Pymol sessions of quality matches

    <Match_Dir>
        Path to original PDB with natively bound target ligand

    <scaffold_file_path>
        Path to file containing list of scaffold PDBs

    <ligand_ID>
        Three-letter ligand code

    --monomer -m
        Monomer matches

    --posfile -p
        Generate new posfiles for iteration
"""
import docopt
import shutil
import os
import re

import pandas as pd
import prody

def copy_matches(pymol_dir, match_dir, final_dir_name='Matches-Final'):
    """
    Copies quality matches from match_dir
    :return:
    """
    # New directory for final match picks
    final_match_dir = os.path.join(os.getcwd(), final_dir_name)
    os.makedirs(final_match_dir, exist_ok=True)

    # Iterate through pymol_dir
    quality_match_list = []
    for pymol_session in os.listdir(pymol_dir):
        if pymol_session.endswith('.pse'):
            quality_match_list.append(pymol_session.replace('-pretty.pse', '.pdb'))
    print('\nFound {0} quality matches in {1}...\n'.format(len(quality_match_list), pymol_dir))

    # Move matches in quality_match_list from match_dir into final_match_dir
    moved_matches = []
    for match in quality_match_list:
        src = os.path.join(match_dir, match)
        if os.path.exists(src):
            dst = os.path.join(final_match_dir, match)
            shutil.copy2(src, dst)
            moved_matches.append(match)

    # Report
    missing_matches = set(quality_match_list) - set(moved_matches)
    if len(missing_matches) == 0:
        print('All quality matches transferred successfully to {0}!'.format(final_match_dir))
    else:
        print('The following matches were not found in {0}:'.format(match_dir))
        for missing in missing_matches:
            print(missing)

    print('\n')

    return final_match_dir


def generate_iteration_inputs(final_match_dir, output_name_prefix='gurobi_constraints-residue_count', monomer=False):
    """
    Generate .csv and .txt for match iteration
    :return:
    """
    df_dict_list = []

    # Generate csv
    for match in os.listdir(final_match_dir):
        if match.endswith('.pdb'):
            conformer_name, constraint_resnums = find_conformer_and_constraint_resnums(match, monomer=monomer)
            scaffold = re.split('_|-|\.', match)[4]
            constraint = '{}-{}'.format(conformer_name, '1_' + '_'.join([str(a) for a in constraint_resnums]))

            temp_dict = {'scaffold': scaffold,
                         'constraint': constraint
                         }

            df_dict_list.append(temp_dict)

    df = pd.DataFrame(df_dict_list)
    df.to_csv('{0}-{1}.csv'.format(output_name_prefix, len(constraint_resnums)), index=False)

    # Get non-redundant list of constraints
    constraint_set = set()
    for index, row in df.iterrows():
        constraint_set.add(row['constraint'])

    # Generate txt
    with open('{0}-{1}.txt'.format(output_name_prefix, len(constraint_resnums)), 'w') as file:
        for constraint in list(constraint_set):
            file.write('{0}\n'.format(constraint))

    # Report
    print('{0} unique constraints found in the current match set\n'.format(len(constraint_set)))
    print('Iteration inputs written to {0}\n'.format(os.getcwd()))


def find_conformer_and_constraint_resnums(pdb_name, monomer=False):
    """
    Generates name of ideal binding motif from match PDB name
    :param pdb_name:
    :return:
    """
    pdb_split = re.split('_|-|\.', pdb_name)

    ligand_name_index = 5
    conformer_id_index = 6

    if len(pdb_split[4]) == 4:
        ligand_name_index += 1
        conformer_id_index += 1

    ligand_name = pdb_split[ligand_name_index]
    conformer_id = pdb_split[conformer_id_index]
    conformer_name = '{}_{}'.format(ligand_name, conformer_id)

    constraint_resnum_block = re.split('-|\.', pdb_name)[1][:-2]
    constraint_resnums = [int(a) for a in constraint_resnum_block.split('_') if a != ''][1:]

    return conformer_name, constraint_resnums


def create_unique_match_ID(match, monomer=False):
    pnc = re.split('_|-|\.', match)
    return '{0}-{1}'.format(pnc[4], pnc[2])


def generate_posfile_things(ligand_ID, scaffold_file, monomer=False):
    """
    Generate new posfiles for match iteration
    :return:
    """

    # Inclulde residues with CA within 8A of ligand
    distance = 9

    # Get relevant positions for all matches
    position_dict = {}

    # New directory for posfiles to be generated
    posfile_path = os.path.join(os.getcwd(), 'posfiles')
    os.makedirs(posfile_path, exist_ok=True)

    # --- Iterate through matches to get match positions --- #
    # Bin posfiles by scaffold and previous match positions. I previously binned on scaffold only, but that blew up the
    # total number of match positions since a single scaffold could potentially accommodate several binding motifs in
    # many different configurations allowed by the inital match posfile.

    for match in os.listdir(final_match_dir):

        if match.endswith('.pdb'):
            unique_mathch_ID = create_unique_match_ID(match, monomer=monomer)
            position_dict[unique_mathch_ID] = set()

    for match in os.listdir(final_match_dir):
        if match.endswith('.pdb'):

            match_prody = prody.parsePDB(os.path.join(final_match_dir, match))
            designable_residues = match_prody.select('(calpha within {0} of resname {1}) and not resname PRO GLY'.format(distance, ligand_ID))
            posfile_list = [atom.getResnum() for atom in designable_residues]

            unique_mathch_ID = create_unique_match_ID(match, monomer=monomer)
            position_dict[unique_mathch_ID] = position_dict[unique_mathch_ID] | set(posfile_list)

    import pprint
    pprint.pprint(position_dict)

    # --- Create posfiles for relevant scaffolds --- #
    scaffold_file = open(scaffold_file)
    scaffold_list = [a.strip() for a in scaffold_file if a]

    relevant_scaffolds = [scaffold for scaffold in scaffold_list if scaffold[:4] in position_dict.keys()]
    for key, value in position_dict.items():
        posfile_name = '{}.pos'.format(key)
        with open(os.path.join(posfile_path, posfile_name), 'w') as posfile:
            posfile.write(' '.join([str(a) for a in value]))


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    pymol_dir = args['<Pymol_Dir>']
    match_dir = args['<Match_Dir>']
    monomer_matches = args['--monomer']
    posfile = args['--posfile']
    ligand_ID = args['<ligand_ID>']
    scaffold_file = args['<scaffold_file_path>']

    # --- Generate match iteration inputs -- #
    final_match_dir = copy_matches(pymol_dir, match_dir)
    generate_iteration_inputs(final_match_dir, monomer=monomer_matches)

    # --- Generate new posfiles for iteration --- #
    if posfile:
        generate_posfile_things(ligand_ID, scaffold_file, monomer_matches)


