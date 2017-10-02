#!/usr/bin/env python3

"""
Filter matches for a given target compound using constraints generated with the BSFF package. 

Usage:
    match_filtering <ligand> <match_PDB_dir> <match_sc_path> <ideal_binding_site_dir> [--monomer] [--csv <csv_path>]

Arguments:     
    <ligand>
        Three letter code of target ligand
        
    <match_PDB_dir>
        Directory containing all matches for a given target compound

    <match_sc_path>
        Path to match_score.sc

    <ideal_binding_site_dir>
        Directory containing ideal binding site PDBs as generated using BSFF 
        under Complete_Matcher_Constraints/Binding_Site_PDBs
        
    <csv_path>
        Path to previously calculated match metrics
        
Options:
    -m --monomer
        Filtering monomer matches
        
    -c --csv
        Import csv at csv_path from previous match set

"""

import docopt
import os
import sys
import re
import pandas as pd
import numpy as np
import prody
import pprint
import shutil
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

    def __init__(self, ligand, match_PDB_dir, ideal_bs_dir, match_sc_path, monomer=False):
        self.ligand = ligand
        self.monomer = monomer
        self.match_PDB_dir = match_PDB_dir
        self.ideal_bs_dir = ideal_bs_dir
        self.ideal_bs_dict = self._import_ideal_binding_sites()
        self.match_score_dict = {line.split()[0]: float(line.split()[1]) for line in open(match_sc_path) if line.split()[0] != 'match_name'}
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
        df_temp = pd.DataFrame(columns=['match_name',
                                        'ligand_shell_eleven',
                                        'interface_CB_contact_percentage',
                                        'motif_shell_CB',
                                        'residue_match_score',
                                        'ligand_match_score',
                                        'min_res_per_chain'
                                        ]
                               )
        return df_temp

    def calculate_CB_stats(self, match_prody, motif_residue_IDs):
        """
        Count CB atoms and all related metrics
        :param match_prody: prody object of matched PDB
        :return: 
        """
        # Calculate number of CB atoms within 11A of ligand
        ligand_shell_eleven = len(match_prody.select('name CB within 11 of resname {}'.format(self.ligand)))
        print('Ligand 11A shell CB count: {}'.format(ligand_shell_eleven))

        # todo: UM_1_Y176W276Q170Q177_1_2BH1_TEP_0001_21-22-25-5_1.pdb oesn't have two chains???
        interface_CB_contact_percentage = 0
        motif_shell_CB = 0

        # Percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are within 6A of any ligand atom
        chains_in_dimer = list(set(match_prody.getChids()) - set('X'))
        if len(chains_in_dimer) == 2:
            interface_cb = match_prody.select('(name CB and chain {}) within 8 of chain {} or\
             (name CB and chain {}) within 8 of chain {}'.format(chains_in_dimer[0], chains_in_dimer[1],
                                                                 chains_in_dimer[1], chains_in_dimer[0]))
            print(len(interface_cb))

            ligand_shell_six = match_prody.select('name CB within 6 of resname {}'.format(self.ligand))
            interface_CB_contact_percentage = len(set(ligand_shell_six.getIndices()) & set(interface_cb.getIndices())) / len(interface_cb)
            print('Interface CB contact percentage: {}'.format(interface_CB_contact_percentage))

        # neighbor bin of motif residues (i.e. number of CB atoms within 8A° of any motif residue CB atom)
        motif_resnums = [res[1] for res in motif_residue_IDs]
        motif_shell_CB = len(match_prody.select('name CB within 8 of (resnum {} or resnum {} or resnum {} or resnum {})'\
                                                .format(motif_resnums[0],
                                                        motif_resnums[1],
                                                        motif_resnums[2],
                                                        motif_resnums[3]
                                                        )
                                                )
                             )

        print('Motif chell CB count: {}'.format(motif_shell_CB))
        return ligand_shell_eleven, interface_CB_contact_percentage, motif_shell_CB

    def calculate_rmsd_stats(self, match_prody, ideal_name, motif_residue_IDs, match_name):
        """
        RMSD things
        :param match_prody: prody of match PDB
        :param ideal_name: name of ideal binding site PDB to be retireved from preloaded dict
        :return: 
        """
        ideal_prody = self.ideal_bs_dict[ideal_name]

        # Calculate RMSD to ideal binding site  (side chains only, not ligand)

        # Atom orders are all messed up in the matched PDBs, need to manually reorder them before aligning coordsets
        ideal_atom_order = ideal_prody.select('resname {}'.format(self.ligand)).getNames()
        match_atom_list = [match_prody.select('resname {} and name {}'.format(self.ligand, atom)) for atom in ideal_atom_order]
        match_atom_coords = np.asarray([atom.getCoords()[0] for atom in match_atom_list])

        # superpose match ligand onto ideal ligand (RMSD should be ~0)
        ideal_ligand = ideal_prody.select('resname {}'.format(self.ligand))
        transformation = prody.calcTransformation(ideal_ligand.getCoords(), match_atom_coords)
        transformed_ideal_prody = prody.applyTransformation(transformation, ideal_prody)

        print('RMSD: {}'.format(prody.calcRMSD(transformed_ideal_prody.select('resname {}'.format(self.ligand)), match_atom_coords)))

        # Debugging
        # prody.writePDB('ideal.pdb', transformed_ideal_prody)
        # prody.writePDB('match.pdb', match_prody.select('resname {}'.format(self.ligand)))

        # Select residues from ideal and get coords
        hv = transformed_ideal_prody.getHierView()
        ideal_residue_prody_list = [res.select('not hydrogen') for res in hv.iterResidues()]

        # Select residues from match and get coords
        motif_residue_prody_list = [match_prody.select('resnum {} and not hydrogen and protein'.format(res_tuple[1])) for res_tuple in motif_residue_IDs]

        # Calculate match score as defined by ME!
        # todo: UM_1_E252W248T285W6_1_3A4U_TEP_0001-11-17-2-22_1.pdb has trouble aligning...
        try:
            residue_match_score = sum([prody.calcRMSD(ideal, match) for ideal, match in zip(ideal_residue_prody_list, motif_residue_prody_list)])
            print('Match score: {}'.format(residue_match_score))
        except:
            residue_match_score = 9999

        # Get ligand match score from matcher_scores.sc
        # todo: sometimes things aren't added to the match_score.sc??
        if match_name.split('.')[0] in list(self.match_score_dict.keys()):
            ligand_match_score = self.match_score_dict[match_name.split('.')[0]]
        else:
            ligand_match_score = 9999

        return residue_match_score, ligand_match_score

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    ligand = args['<ligand>']
    match_PDB_dir = args['<match_PDB_dir>']
    ideal_bs_dir = args['<ideal_binding_site_dir>']
    monomer = args['--monomer']
    match_sc_path = args['<match_sc_path>']

    res_one_to_three = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
                        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
                        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
                        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}

    if args['--csv']:
        df = pd.read_csv(args['<csv_path>'])

    else:
        filter = Filter_Matches(ligand, match_PDB_dir, ideal_bs_dir, match_sc_path, monomer=monomer)
        pprint.pprint(filter.ideal_bs_dict)
        row_list = []

        # Calculate stats for each matched PDB
        for matched_PDB in pdb_check(match_PDB_dir):
            match_name = os.path.basename(os.path.normpath(matched_PDB))
            print(match_name)
            match_prody = prody.parsePDB(matched_PDB)

            # Parse matched_PDB to get ideal binding site name and residues
            pnc = re.split('_|-|\.', os.path.basename(os.path.normpath(matched_PDB)))

            # Ideal Binding Site Name
            ideal_binding_site_name = '{}_{}-{}-{}-{}-{}'.format(pnc[5], pnc[6], pnc[7], pnc[8], pnc[9], pnc[10])
            print(ideal_binding_site_name)
            # Motif Residues in Matched PDB
            motif_residue_ID_list = [a for a in re.split('(\D+)', pnc[2]) if a != '']
            motif_residue_IDs = [(res_one_to_three[motif_residue_ID_list[indx]], motif_residue_ID_list[indx + 1]) for indx in range(0, len(motif_residue_ID_list), 2)]

            # Calculate number of CB atoms within 10A of ligand
            # Percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are within 6A of any ligand atom

            ligand_shell_eleven, interface_CB_contact_percentage, motif_shell_CB = filter.calculate_CB_stats(match_prody, motif_residue_IDs)

            # Calculate match score as defined by Roland
            # Calculate RMSD to ideal binding site  (side chains only, not ligand)

            residue_match_score, ligand_match_score = filter.calculate_rmsd_stats(match_prody, ideal_binding_site_name, motif_residue_IDs, match_name)

            # minimum number of motif residues per chain
            # todo: accomodate cases where all residues are on one chain! Currently returns 4 b/c list of res.getChIDs() is used to determine this
            motif_resnums = [res[1] for res in motif_residue_IDs]
            motif_residues = [match_prody.select('resnum {}'.format(motif_resnum)) for motif_resnum in motif_resnums]
            motif_residue_chain_list = [res.getChids()[0] for res in motif_residues]

            # todo: UM_1_D267F289Y271Q279_1_2BH1_TEP_0001-10-18-21-25_1 is empty??
            try:
                min_res_per_chain = min([motif_residue_chain_list.count(chain) for chain in (set(motif_residue_chain_list) - set('X'))])
                if min_res_per_chain == 4:
                    min_res_per_chain = 0
            except:
                min_res_per_chain = -1

            print(min_res_per_chain)
            print('\n')

            # Aggragate results
            row_dict = {'match_name': match_name,
                        'ligand_shell_eleven': ligand_shell_eleven,
                        'interface_CB_contact_percentage': interface_CB_contact_percentage,
                        'motif_shell_CB': motif_shell_CB,
                        'residue_match_score': residue_match_score,
                        'ligand_match_score': ligand_match_score,
                        'min_res_per_chain': min_res_per_chain
                        }
            row_list.append(row_dict)

        # Return passing matcher results
        df = pd.DataFrame(row_list)
        df.set_index(['match_name'], inplace=True)
        df.to_csv('RESULTS.csv')

    # Let's say take top 5% of hits, for each metric, passing matcher results have to be in all top 5%
    df.set_index(['match_name'], inplace=True)
    percent_cutoff = 0.2
    cut_index = int(len(df) * percent_cutoff)

    # Ascending = True or False
    ascending_dict = {'ligand_match_score': True,
                      'ligand_shell_eleven': False,
                      'interface_CB_contact_percentage': False,
                      'motif_shell_CB': False,
                      'residue_match_score': True,
                      'min_res_per_chain': True
                      }

    set_list = []
    pprint.pprint(list(ascending_dict.keys()))
    for index, column in df.iteritems():
        if index in list(ascending_dict.keys()):
            print(index)
            # print(column.sort_values(ascending=ascending_dict[index])[:cut_index])
            set_list.append(set(column.sort_values(ascending=ascending_dict[index]).index.tolist()[:cut_index]))

    # Get all matches with min_res_per_chain = 0
    min_res_none_set = set(df.groupby(['min_res_per_chain']).get_group(0).index.tolist())

    initial_set = set.intersection(*set_list)
    final_set = initial_set - min_res_none_set

    pprint.pprint(final_set)
    print(len(final_set))

    # Move matches that pass filters into a new directory
    filtered_matches_dir = 'Matches-Filtered'
    os.makedirs(filtered_matches_dir, exist_ok=True)
    for match in list(final_set):
        shutil.copy(os.path.join(match_PDB_dir, match), filtered_matches_dir)