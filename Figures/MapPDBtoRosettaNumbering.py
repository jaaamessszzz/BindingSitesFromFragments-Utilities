#!/usr/bin/env python3

"""
Create pickles that map Rosetta numbering in clean scaffold to PDB numbering
Walk through project directories and only generate pickles for protein chains in clean PDB used for applications.

Usage:
    MapPDBtoRosettaNumbering map <compound_dir>
    MapPDBtoRosettaNumbering report <map_dir> <designable_pickle>

Arguments:
    <compound_dir>      Top level directory containing BSFF project directories
    <fasta_dir>         Directory containing FASTA files for designable positions in scaffolds
"""
import os
import sys
import shutil
import pickle
from pprint import pprint

import prody
import docopt

def clean_pdb_mapping(input_pdb, param_pdb, ligand_code=None, protein_chains=None, ligand_ids=[], hetatm_ids=None):
    """
    Clean PDB. Roland: MSE to MET records, CSE to CYS records, discarding alternative conformations, setting all atom
    occupancies to 1.0, discarding all residues with any missing main chain atom(s), removing all ligands but the
    target, and ensuring internal Rosetta numbering

    NOTE: This is a version of the function from ligand_interface_energy.py modified to keep track of clean->PDB

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

    # Rosetta to PDB numbering
    rosetta_to_pdb_mapping = {}

    for chain in hv:
        for residue in chain:
            if (residue.getChid(), residue.getResnum()) in ligand_ids:
                param_pdb_prody = prody.parsePDB(param_pdb)

            # Check Backbone atoms, else don't add to output atomgroup
            elif all(atom in residue.getNames() for atom in ['N', 'C', 'CA']) and residue.getResname() in canon_residues:
                rosetta_to_pdb_mapping[res_count] = f'{residue.getResnum()}{residue.getChid()}'
                residue.setResnum(res_count)
                res_count += 1

                if residue.getResname() in ['MSE', 'SEC']:
                    residue = _fix_mse_sec(residue, residue.getResname())

                residue.setOccupancies([1] * len(residue))
                output_atoms = _add_to_output(output_atoms, residue)

            elif residue.getResname() in ['A', 'T', 'C', 'G', 'U']:
                rosetta_to_pdb_mapping[res_count] = f'{residue.getResnum()}{residue.getChid()}'
                residue.setResnum(res_count)
                res_count += 1
                output_atoms = _add_to_output(output_atoms, residue)

            elif residue.getResname() in ['CA', 'SO4']:
                rosetta_to_pdb_mapping[res_count] = f'{residue.getResnum()}{residue.getChid()}'
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

    return output_atoms, removed_residues, rosetta_to_pdb_mapping


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    if args['<compound_dir>']:
        compound_dir = args['<compound_dir>']
        dump_dir = os.path.join(os.getcwd(), 'RosettaResnumMapping')
        os.makedirs(dump_dir, exist_ok=True)

        for ccid in os.listdir(compound_dir):
            project_root = os.path.join(compound_dir, ccid)
            if not ccid.startswith('_') and os.path.isdir(project_root):
                # All relevant information is in DesignInputs
                design_inputs_path = os.path.join(project_root, 'DesignInputs')
                design_inputs_list = os.listdir(design_inputs_path)

                complex_pdb_gz = [a for a in design_inputs_list if a.endswith('.pdb.gz')][0]
                clean_pdb = [a for a in design_inputs_list if a.endswith('-clean.pdb')][0]
                param_pdb = [a for a in design_inputs_list if a.endswith('_0001.pdb')][0]

                pdbid = complex_pdb_gz.split('.')[0]
                complex_pdb_prody = prody.parsePDB(os.path.join(design_inputs_path, complex_pdb_gz))

                # Get Chains from clean PDB
                clean_pdb_prody = prody.parsePDB(os.path.join(design_inputs_path, clean_pdb))
                clean_chains = [chain.getChid() for chain in clean_pdb_prody.getHierView()]

                # Generate scaffold-specific pickle with same name as clean PDB
                output_atoms, removed_residues, rosetta_to_pdb_mapping = clean_pdb_mapping(complex_pdb_prody, param_pdb, protein_chains=clean_chains)
                pprint(rosetta_to_pdb_mapping)
                with open(os.path.join(dump_dir, f'{clean_pdb[:4]}.pickle'), 'wb') as output:
                    pickle.dump(rosetta_to_pdb_mapping, output)

    if args['<map_dir>']:
        design_positions = args['<designable_pickle>']
        resnum_mapping_dir = args['<map_dir>']

        design_positions_dict = pickle.load(open(design_positions, 'rb'))
        pprint(design_positions_dict)

        for file in os.listdir(resnum_mapping_dir):
            if file.endswith('.pickle'):
                complex_name = file.split('.')[0]
                complex_mapping = pickle.load(open(os.path.join(resnum_mapping_dir, file), 'rb'))
                complex_designable = design_positions_dict[complex_name]['positions']
                pdb_numbering = [complex_mapping[position] for position in complex_designable]

                print(complex_name)
                print(' ,'.join(pdb_numbering))
