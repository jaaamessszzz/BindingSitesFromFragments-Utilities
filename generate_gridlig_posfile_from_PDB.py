#!/usr/bin/env python3
"""
Generates inputs for the matcher given an input pdb containing a WT binding site with ligand bound.

The purpose of this script is to easily generate the inputs required for the matcher to test whether we can recapitulate
native binding sites for arbitrary small molecules with the BSFF protocol.

Input...
* PDB with ligand
* Three letter code of ligand (passed through cmd line)
* Resnum of ligand (passed through cmd line)

Output...
* Cleaned, renumbered PDB
* Gridlig file (Roland: unit size of 0.1A with a 10A buffer on each side)
* Posfile file (Roland: all residues (in Rosetta numbering) that are within 15A° distance of the given other chain)

I'm anticipating that we will be mostly dealing with monomers, so Posfile residues wil be determined to be residues 
within 15A of the natively bound ligand.

Usage:
    generate_gridlig_posfile_from_PDB <input_PDB> <ligand> <resnum>

Arguments:    
    <input_PDB>
        Path to original PDB with natively bound target ligand
        
    <ligand>
        Three letter code of target ligand
        
    <resnum>
        Residue number of target number within <input_PDB>
"""
import docopt
import prody
import sys
import os
import numpy as np

def clean_pdb(input_pdb, ligand_code, resnum):
    """
    Clean PDB. Roland: MSE to MET records, CSE to CYS records, discarding alternative conformations, setting all atom 
    occupancies to 1.0, discarding all residues with any missing main chain atom(s), removing all ligands but the 
    target, and ensuring internal Rosetta numbering
    :return: 
    """
    # Apparently this selects Altloc A if there are alternate locations at all... makes things easy
    cleanish_pdb = input_pdb.select('(protein or hetero) and not water').copy()
    hv = cleanish_pdb.getHierView()

    # Output atomgroup
    output_atoms = prody.atomgroup.AtomGroup('Output')

    # Residue count
    res_count = 1

    for chain in hv:
        for residue in chain:
            if residue.getResnum() == resnum and residue.getResname() == ligand_code:

                # Set Rosetta Numbering
                residue.setResnum(res_count)
                res_count += 1

                # Add to output atomgroup
                print(residue)
                output_atoms = _add_to_output(output_atoms, residue)

            # Check Backbone atoms, else don't add to output atomgroup
            elif all(atom in residue.getNames() for atom in ['N', 'C', 'CA']):

                # Set Rosetta Numbering
                residue.setResnum(res_count)
                res_count += 1

                # Check and correct for MSE/SEC
                if residue.getResname() in ['MSE', 'SEC']:
                    residue = _fix_mse_sec(residue, residue.getResname())

                # Set all occupancies to 1
                residue.setOccupancies([1] * len(residue))

                # Add to output atomgroup
                print(residue)
                output_atoms = _add_to_output(output_atoms, residue)

            # Tossed HETATM or residues
            else:
                print('Removed {}'.format(residue))

    return output_atoms

def generate_posfile_info(cleaned_pdb, ligand_code):
    """
    Roland: all residues (in Rosetta numbering) that are within 15A° distance of the given other chain
    :param cleaned_pdb: cleaned_pdb from clean_pdb()
    :param ligand_code: three letter ligand code
    :return: list of designable resnums in cleaned pdb
    """
    designable_residues = cleaned_pdb.select('calpha within 15 of resname {}'.format(ligand_code))
    return [atom.getResnum() for atom in designable_residues]

def generate_gridlig_info(cleaned_pdb, ligand_code):
    """
    Roland: unit size of 0.1A with a 10A buffer on each side
    :param cleaned_pdb: cleaned_pdb from clean_pdb()
    :param ligand_code: three letter ligand code
    :return: lower corner of gridlig
    """
    ligand = cleaned_pdb.select('resname {}'.format(ligand_code))
    ligand_center = prody.calcCenter(ligand)
    return (ligand_center - np.array([5,5,5])).tolist()


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
    seleno_index = [e for e in res_elements].index('SE')

    # Set SE to S
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

def main():
    args = docopt.docopt(__doc__)
    input_pdb = prody.parsePDB(args['<input_PDB>'])
    ligand_code = args['<ligand>'].upper()
    resnum = int(args['<resnum>'])

    print(input_pdb, ligand_code, resnum)

    # Generate cleaned PDB
    cleaned_pdb = clean_pdb(input_pdb, ligand_code, resnum)
    # todo: get rid of the fucking remarks in the output PDB
    prody.writePDB('{}_clean.pdb.gz'.format(args['<input_PDB>'].upper()), cleaned_pdb)
    scaffold_only = cleaned_pdb.select('protein')
    prody.writePDB('{}_scaffold.pdb.gz'.format(args['<input_PDB>'].upper()), scaffold_only)

    # Generate Posfile
    posfile_list = generate_posfile_info(cleaned_pdb, ligand_code)
    with open('{}.pdb.pos'.format(args['<input_PDB>'].upper()), 'w') as posfile:
        posfile.write(' '.join([str(a) for a in posfile_list]))

    # Generate Gridlig... seems silly to do if I'm not specifying allowed voxels...
    gridlid_lower_corner = generate_gridlig_info(cleaned_pdb, ligand_code)
    with open('{}.pdb.gridlig'.format(args['<input_PDB>'].upper()), 'w') as gridlig:
        gridlig.write('NAME: gridlig\n')
        gridlig.write('BASE: {}\n'.format(' '.join([str(coord) for coord in gridlid_lower_corner])))
        gridlig.write('SIZE: 100 100 100\n')
        gridlig.write('LENGTH: 0.1 0.1 0.1')

if __name__ == "__main__":
    main()