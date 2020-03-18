#! /usr/bin/python3

"""
Rank residues at protein-ligand interface by interaction energy

Takes a cleaned PDB in Rosetta numbering and outputs a dataframe for interface residues ranked by REU

Usage:
  ligand_interface_energy <pdb_path> <ligand_code> <params_file> <param_pdb> <protein_chains> (<ligand_chain> <ligand_position>)... [options]

Arguments:
  <pdb_path>            Path to PDB for scoring
  <ligand_code>         Chemical component identifier for ligand
  <params_file>         Path to Rosetta params file for ligand
  <param_pdb>           Path to PDB created by molfile_to_params
  <protein_chains>      Chain IDs for protein chains to keep
  <ligand_chain>        Chain ID for ligand
  <ligand_position>     Position of ligand to calculate interface energies

Options:
  --clean, -c                               Clean input PDB before assessing energies
"""


import docopt
import io
import os
import pandas as pd
import prody
from pprint import pprint
from rdkit import Chem

import pyrosetta
from pyrosetta import rosetta

# todo: make and install a package with all these utility functions...
def clean_pdb(input_pdb, param_pdb, ligand_code=None, protein_chains=None, ligand_ids=None, hetatm_ids=None):
    """
    Clean PDB. Roland: MSE to MET records, CSE to CYS records, discarding alternative conformations, setting all atom
    occupancies to 1.0, discarding all residues with any missing main chain atom(s), removing all ligands but the
    target, and ensuring internal Rosetta numbering

    :param input_pdb: prody of PDB to clean
    :param ligand_code: chemical component identifier for ligand to keep
    :param ligand_positions: position of ligand in input_pdb to keep
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
    cleanish_pdb = input_pdb.select(f'(protein or hetero or nucleic) and not water and chain {" ".join([chain for chain in protein_chains])}').copy()
    hv = cleanish_pdb.getHierView()

    output_atoms = prody.atomgroup.AtomGroup('Output')
    param_pdb_prody = None
    res_count = 1
    removed_residues = []

    for chain in hv:
        for residue in chain:
            if (residue.getChid(), residue.getResnum()) in ligand_ids:
                param_pdb_prody = prody.parsePDB(param_pdb)

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

            elif residue.getResname() in ['CA', 'SO4']:
                residue.setResnum(res_count)
                res_count += 1
                residue.setOccupancies([1] * len(residue))
                output_atoms = _add_to_output(output_atoms, residue)

            else:
                print('Removed {}'.format(residue))
                removed_residues.append('{0}{1}'.format(residue.getResname(), residue.getResnum()))

    # Append ligand to end of chain
    if param_pdb_prody:
        param_pdb_prody.setResnums([res_count] * len(param_pdb_prody))
        output_atoms = _add_to_output(output_atoms, param_pdb_prody)

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
    bsff_sfxn.set_weight(rosetta.core.scoring.hbond_bb_sc, 1)

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
    pprint(args)

    input_pdb_path = args['<pdb_path>']
    ligand_cci = args['<ligand_code>']
    ligand_params = args['<params_file>']
    param_pdb = args['<param_pdb>']
    protein_chains = args['<protein_chains>']

    ligand_chains = [a for a in args['<ligand_chain>']] if type(args['<ligand_chain>']) is list else list(args['<ligand_chain>'])
    ligand_positions = [int(a) for a in args['<ligand_position>']] if type(args['<ligand_position>']) is list else list(args['<ligand_position>'])
    ligand_ids = [(chain, pos) for chain, pos in zip(ligand_chains, ligand_positions)]
    pyrosetta.init(f'-extra_res_fa {ligand_params}')

    if args['--clean']:
        prody_intermediate = prody.parsePDB(input_pdb_path)
        cleaned_prody_intermediate, removed_residues = clean_pdb(prody_intermediate, param_pdb, ligand_code=ligand_cci, protein_chains=protein_chains, ligand_ids=ligand_ids)


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
        ligand_position = cleaned_prody_intermediate.select(f'resname {ligand_cci} and chain X').getResnums()[0] # Jank
        print(f'New ligand position is {ligand_cci}:{ligand_position}')

        # Generate params file for ligand


    else:
        complex_pose = rosetta.core.import_pose.pose_from_file(input_pdb_path)

    # Relax
    sfxn = rosetta.core.scoring.get_score_function()
    fast_relax = rosetta.protocols.relax.FastRelax(sfxn, 5, 'MonomerRelax2019')
    fast_relax.constrain_relax_to_native_coords(True)
    fast_relax.apply(complex_pose)
    complex_pose.dump_pdb(f'{os.path.basename(input_pdb_path)[:4]}-relaxed.pdb')

    interface_df = calculate_ligand_interfaceE(complex_pose, ligand_cci, ligand_position, ligand_params)
    interface_df.sort_values('interE', inplace=True)
    interface_df.to_csv(f'{os.path.basename(input_pdb_path)[:4]}-LigandInterfaceEnergy.csv')
    complex_pose.dump_pdb(f'{os.path.basename(input_pdb_path)[:4]}-scored.pdb')
