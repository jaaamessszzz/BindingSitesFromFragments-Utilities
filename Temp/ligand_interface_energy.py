#! /usr/bin/python3

"""
Rank residues at protein-ligand interface by interaction energy

Takes a cleaned PDB in Rosetta numbering and outputs a dataframe for interface residues ranked by REU

Usage:
  ligand_interface_energy <pdb_path> <ligand_code> <ligand_position> <params_file> [options]

Arguments:
  <pdb_path>            Path to PDB for scoring
  <ligand_code>         Chemical component identifier for ligand
  <ligand_position>     Position of ligand to calculate interface energies
  <params_file>         Path to Rosetta params file for ligand

Options:
  --clean, -c           Clean input PDB before assessing energies

"""


import docopt
import io
import os
import pandas as pd
import prody

import pyrosetta
from pyrosetta import rosetta

# todo: make and install a package with all these utility functions...
def clean_pdb(input_pdb, ligand_code=None, resnum=None):
    """
    Clean PDB. Roland: MSE to MET records, CSE to CYS records, discarding alternative conformations, setting all atom
    occupancies to 1.0, discarding all residues with any missing main chain atom(s), removing all ligands but the
    target, and ensuring internal Rosetta numbering

    :param input_pdb: prody of PDB to clean
    :param ligand_code: chemical component identifier for ligand to keep
    :param resnum: position of ligand in input_pdb to keep
    :return:
    """

    canon_residues = ['ALA', 'CYS', 'SEC', 'ASP', 'GLU',
                      'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                      'LEU', 'MET', 'MSE', 'ASN', 'PRO',
                      'GLN', 'ARG', 'SER', 'THR', 'VAL',
                      'TRP', 'TYR']

    # Necessary due to prody things...
    def _add_to_output(output_atoms, residue):
        if len(output_atoms) != 0:
            output_atoms = output_atoms + residue.copy()
        else:
            output_atoms = residue.copy()
        return output_atoms

    # Pirated from Generate_Motif_Residues
    def _fix_mse_sec(representative_residue, resname):
        """
        Takes MSE/SEC and turns it into MET/CYS
        :param representative_residue: representative_residue with Se
        :return: representative_residue with S
        """
        # Find index of SE
        res_elements = representative_residue.getElements()

        # If SE is missing due to missing atoms, just return the residue as is
        if 'SE' not in res_elements:
            return representative_residue

        # Set SE to S
        seleno_index = [e for e in res_elements].index('SE')
        res_elements[seleno_index] = 'S'

        # Set elements to MET
        representative_residue.setElements(res_elements)

        # Set resnames MSE>MET, SEC>CYS
        # Set seleno atom name to met/cys atom sulfur atom name
        if resname == 'MSE':
            representative_residue.setResnames(['MET'] * len(representative_residue))
            res_atom_names = representative_residue.getNames()
            res_atom_names[seleno_index] = 'SD'
            representative_residue.setNames(res_atom_names)
        elif resname == 'SEC':
            representative_residue.setResnames(['CYS'] * len(representative_residue))
            res_atom_names = representative_residue.getNames()
            res_atom_names[seleno_index] = 'SG'
            representative_residue.setNames(res_atom_names)
        else:
            raise Exception('Either MSE or SEC had to be set!')

        return representative_residue

    # This selects Altloc A if there are alternate locations at all... makes things easy
    cleanish_pdb = input_pdb.select('(protein or hetero or nucleic) and not water').copy()
    hv = cleanish_pdb.getHierView()

    output_atoms = prody.atomgroup.AtomGroup('Output')
    res_count = 1
    removed_residues = []

    for chain in hv:
        for residue in chain:
            if residue.getResnum() == resnum and residue.getResname() == ligand_code:
                residue.setResnum(res_count)
                res_count += 1
                output_atoms = _add_to_output(output_atoms, residue)

            # Check Backbone atoms, else don't add to output atomgroup
            elif all(atom in residue.getNames() for atom in ['N', 'C', 'CA']) and residue.getResname() in canon_residues:
                residue.setResnum(res_count)
                res_count += 1

                if residue.getResname() in ['MSE', 'SEC']:
                    residue = _fix_mse_sec(residue, residue.getResname())

                residue.setOccupancies([1] * len(residue))
                output_atoms = _add_to_output(output_atoms, residue)

            elif residue.getResname() in ['A', 'T', 'C', 'G', 'U']:
                residue.setResnum(res_count)
                res_count += 1
                output_atoms = _add_to_output(output_atoms, residue)

            else:
                print('Removed {}'.format(residue))
                removed_residues.append('{0}{1}'.format(residue.getResname(), residue.getResnum()))

    return output_atoms, removed_residues

def calculate_ligand_interfaceE(input_pdb, ligand_cci, ligand_position, ligand_params):
    """
    Iterate over residues in a pose and return a sorted dataframe of interface energies

    :param input_pdb: pyrosetta pose of cleaned PDB
    :param ligand_cci: chemical component identifier of target ligand
    :param ligand_position: position of
    :param ligand_params:
    :return:
    """

    # Infrastructure
    score_dict_list = list()

    # Custom score function using select 2b energies
    bsff_sfxn = rosetta.core.scoring.ScoreFunction()
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_atr, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_rep, 0.55)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_sol, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_elec, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.hbond_sc, 1)

    bsff_sfxn(input_pdb)
    sfxn_weights = bsff_sfxn.weights()
    edges = input_pdb.energies().energy_graph()

    for position in range(1, input_pdb.size() + 1):
        if position == ligand_position:
            continue

        current_edge = edges.find_energy_edge(position, ligand_position)
        interface_reu = current_edge.dot(sfxn_weights) if current_edge is not None else 0

        score_dict = {'position': position,
                      'resname': input_pdb.residue(position).name3(),
                      'interE': interface_reu}

        score_dict_list.append(score_dict)

    return pd.DataFrame(score_dict_list)


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    input_pdb_path = args['<pdb_path>']
    ligand_cci = args['<ligand_code>']
    ligand_position = int(args['<ligand_position>'])
    ligand_params = args['<params_file>']

    pyrosetta.init(f'-extra_res_fa {ligand_params}')

    if args['--clean']:
        prody_intermediate = prody.parsePDB(input_pdb_path)
        cleaned_prody_intermediate, removed_residues = clean_pdb(prody_intermediate, ligand_code=ligand_cci, resnum=ligand_position)
        prody.writePDB(f'{os.path.basename(input_pdb_path)[:4]}-clean.pdb', cleaned_prody_intermediate)
        prody_stream = io.StringIO()
        prody.writePDBStream(prody_stream, cleaned_prody_intermediate)
        complex_pose = rosetta.core.pose.Pose()
        prody_stream.seek(0)
        rosetta.core.import_pose.pose_from_pdbstring(complex_pose, prody_stream.read())

        # todo: get position of ligand in new clean residue numbering
        print(ligand_cci)
        print(cleaned_prody_intermediate.select(f'resname {ligand_cci}'))
        print(set(cleaned_prody_intermediate.getResnames()))
        ligand_position = cleaned_prody_intermediate.select(f'resname {ligand_cci}').getResnums()[0] # Jank
        print(f'New ligand position is {ligand_cci}:{ligand_position}')

    else:
        complex_pose = rosetta.core.import_pose.pose_from_file(input_pdb_path)

    interface_df = calculate_ligand_interfaceE(complex_pose, ligand_cci, ligand_position, ligand_params)
    interface_df.sort_values('interE', inplace=True)
    interface_df.to_csv(f'{os.path.basename(input_pdb_path)[:4]}-LigandInterfaceEnergy.csv')