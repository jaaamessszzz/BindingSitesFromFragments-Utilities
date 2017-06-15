#!/usr/bin/env python3

"""
Filter matches for a given target compound using constraints generated with the 
BSFF package. 

Usage:
    match_filtering <ligand> <match_PDB_dir> <ideal_binding_site_dir> [--monomer]

Arguments:     
    <ligand>
        Three letter code of target ligand
        
    <match_PDB_dir>
        Directory containing all matches for a given target compound
    
    <ideal_binding_site_dir>
        Directory containing ideal binding site PDBs as generated using BSFF 
        under Complete_Matcher_Constraints/Binding_Site_PDBs
        
Options:
    -m --monomer
        Filtering monomer matches

"""

import docopt
import os
import sys
import re
import pandas as pd
import numpy as np
import prody
import pprint
from utils import pdb_check

class Filter_Matches:
    """
    Everything required to filter matches for a given target compound using
    constraints generated with the BSFF package.
    
    Adapted from Roland's biosensor design protocol:
    * CB within shell around ligand (10A)
    * percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are
    within 6A° of any ligand atom
    * RMSD to defined motif
    * match score (i.e. sum of RMSDs between all pairs of ligand placements from individual motif
    residues)
    * neighbor bin of motif residues (i.e. number of C" atoms within 8A° of any motif residue C"
    atom)
    * minimum number of motif residues per chain
    """

    def __init__(self, ligand, match_PDB_dir, ideal_bs_dir, monomer=False):
        self.ligand = ligand
        self.monomer = monomer
        self.match_PDB_dir = match_PDB_dir
        self.ideal_bs_dir = ideal_bs_dir
        self.ideal_bs_dict = self._import_ideal_binding_sites()
        self.df = self._set_up_dataframe()

    def _import_ideal_binding_sites(self):
        """
        Import all ideal binding site PDBs whose constraints were used to find
        successful matches.
        :return: 
        """
        ideal_bs_prody_dict = {os.path.basename(os.path.normpath(binding_site_PDB)).split('.')[0]: prody.parsePDB(binding_site_PDB)
                               for binding_site_PDB in pdb_check(self.ideal_bs_dir)}

        return ideal_bs_prody_dict

    def _set_up_dataframe(self):
        """
        Set up a Pandas Dataframe to record all relevant match metrics
        :return: 
        """
        df_temp = pd.DataFrame()

        return df_temp

    def calculate_CB_stats(self, match_prody):
        """
        Count CB atoms and all related metrics
        :param match_prody: prody object of matched PDB
        :return: 
        """
        # Calculate number of CB atoms within 10A of ligand
        ligand_shell_ten = match_prody.select('name CB within 10 of resname {}'.format(self.ligand))
        print(len(ligand_shell_ten))

        # Percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are within 6A of any ligand atom
        chains_in_dimer = list(set(match_prody.getChids()) - set('X'))
        if len(chains_in_dimer) == 2:
            interface_cb = match_prody.select('(name CB and chain {}) within 8 of chain {} or\
             (name CB and chain {}) within 8 of chain {}'.format(chains_in_dimer[0], chains_in_dimer[1],
                                                                 chains_in_dimer[1], chains_in_dimer[0]))
            print(len(interface_cb))

            ligand_shell_six = match_prody.select('name CB within 6 of resname {}'.format(self.ligand))
            interface_CB_contact_percentage = len(set(ligand_shell_six.getIndices()) & set(interface_cb.getIndices())) / len(interface_cb)
            print(interface_CB_contact_percentage)

    def calculate_rmsd_stats(self, match_prody, ideal_name, motif_residue_IDs):
        """
        RMSD things
        :param match_prody: prody of match PDB
        :param ideal_name: name of ideal binding site PDB to be retireved from preloaded dict
        :return: 
        """
        ideal_prody = self.ideal_bs_dict[ideal_name]
        print(ideal_prody)
        print(match_prody)

        # Calculate RMSD to ideal binding site  (side chains only, not ligand)

        # Atom orders are all messed up in the matched PDBs, need to manually reorder them before aligning coordsets
        ideal_atom_order = ideal_prody.select('resname {}'.format(self.ligand)).getNames()
        match_atom_list = [match_prody.select('resname {} and name {}'.format(self.ligand, atom)) for atom in ideal_atom_order]
        match_atom_coords = np.asarray([atom.getCoords()[0] for atom in match_atom_list])

        # superpose match ligand onto ideal ligand (RMSD should be ~0)
        ideal_ligand = ideal_prody.select('resname {}'.format(self.ligand))
        transformation = prody.calcTransformation(ideal_ligand.getCoords(), match_atom_coords)
        transformed_ideal_prody = prody.applyTransformation(transformation, ideal_prody)

        print('RMSD')
        print(prody.calcRMSD(transformed_ideal_prody.select('resname {}'.format(self.ligand)), match_atom_coords))

        # Debugging
        prody.writePDB('ideal.pdb', transformed_ideal_prody)
        prody.writePDB('match.pdb', match_prody.select('resname {}'.format(self.ligand)))

        # Select residues from ideal and get coords
        print([res for res in transformed_ideal_prody])
        ideal_res_coord_list = []

        # Select residues from match and get coords
        print(motif_residue_IDs)
        match_res_coord_list = []

        # calc rmsd between superposed match and ideal

        # Calculate match score as defined by Roland

        # sum(RMSD For each residue in superposed match and ideal)



if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    ligand = args['<ligand>']
    match_PDB_dir = args['<match_PDB_dir>']
    ideal_bs_dir = args['<ideal_binding_site_dir>']
    monomer = args['--monomer']

    res_one_to_three = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
                        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
                        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
                        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}

    filter = Filter_Matches(ligand, match_PDB_dir, ideal_bs_dir, monomer=monomer)

    pprint.pprint(filter.ideal_bs_dict)

    # Calculate stats for each matched PDB
    for matched_PDB in pdb_check(match_PDB_dir):
        print(os.path.basename(os.path.normpath(matched_PDB)))
        match_prody = prody.parsePDB(matched_PDB)

        # Calculate number of CB atoms within 10A of ligand
        # Percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are within 6A of any ligand atom

        filter.calculate_CB_stats(match_prody)

        # Calculate match score as defined by Roland
        # Calculate RMSD to ideal binding site  (side chains only, not ligand)

        # Parse matched_PDB to get ideal binding site name and residues
        pnc = re.split('_|-|\.', os.path.basename(os.path.normpath(matched_PDB)))

        # Ideal Binding Site Name
        ideal_binding_site_name = '{}_{}-{}-{}-{}-{}'.format(pnc[5], pnc[6], pnc[7], pnc[8], pnc[9], pnc[10])
        print(ideal_binding_site_name)
        # Motif Residues in Matched PDB
        motif_residue_ID_list = [a for a in re.split('(\D+)', pnc[2]) if a != '']
        motif_residue_IDs = [(res_one_to_three[motif_residue_ID_list[indx]], motif_residue_ID_list[indx + 1]) for indx in range(0, len(motif_residue_ID_list), 2)]

        filter.calculate_rmsd_stats(match_prody, ideal_binding_site_name, motif_residue_IDs)

        # neighbor bin of motif residues (i.e. number of CB atoms within 8A° of any motif residue CB atom)
        # minimum number of motif residues per chain

        break
    # Aggragate results
    # Return passing matcher results
    # Let's say take top 5% of hits, for each metric, passing matcher results have to be in all top 5%