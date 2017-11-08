#!/usr/bin/env python3

"""
Filter matches for a given target compound using constraints generated with the BSFF package. 

Usage:
    match_filtering <ligand> <match_PDB_dir> <ideal_binding_site_dir> <gurobi_solutions_dir> <match_sc_path> [--monomer] [--csv <csv_path>]

Arguments:     
    <ligand>
        Three letter code of target ligand
        
    <match_PDB_dir>
        Directory containing all matches for a given target compound

    <match_sc_path>
        Path to match_score.sc
    
    <gurobi_solutions_dir>
        Path to directory containing csvs with Gurobi solutions and associated motif scores
    
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
from ast import literal_eval
from utils import pdb_check, minimum_contact_distance
from pathos.multiprocessing import ProcessingPool as Pool

# todo: filter out PDBs that contains gigantic AF chains
# todo: add support for iterative filtering, we'll really only care about match_score, ligand_score, binding_motif_energy

class Filter_Matches:
    """
    Everything required to filter matches for a given target compound using
    constraints generated with the BSFF package.
    
    Adapted from Roland's biosensor design protocol:
    * CB within shell around ligand (10A)
    * percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are within 6A° of any ligand atom
    * RMSD to defined motif
    * match score (i.e. sum of RMSDs between all pairs of ligand placements from individual motif residues)
    * neighbor bin of motif residues (i.e. number of C" atoms within 8A° of any motif residue C" atom)
    * minimum number of motif residues per chain
    """

    def __init__(self, ligand, match_PDB_dir, ideal_bs_dir, match_sc_path, monomer=False):
        self.ligand = ligand
        self.monomer = monomer
        self.match_PDB_dir = match_PDB_dir
        self.ideal_bs_dir = ideal_bs_dir
        self.ideal_bs_dict = self._import_ideal_binding_sites()
        self.match_score_dict = self.setup_match_score_dict(match_sc_path)
        self.df = self._set_up_dataframe()

    def setup_match_score_dict(self, match_sc_path):
        """
        Convert match_scores.sc into a dict with "Match": "score" key value pairs. This method is necessary since all
        match scores are written to match_scores.sc even if output_matches_per_group == 1
        :param match_sc_path: path to match_scores.sc
        :return: 
        """
        match_score_dict_list = []
        with open(match_sc_path, 'r') as match_sc_contents:
            for line in match_sc_contents:
                if line.strip() not in [None, '', '\n']:
                    print('line:', line)
                    if line.split()[0].startswith('UM_'):
                        match_score_dict_list.append({'match_name': line.split()[0], 'match_score': float(line.split()[1])})

        df = pd.DataFrame(match_score_dict_list, columns=['match_name', 'match_score'])
        match_score_dict = {}

        match_PDBs = [match.split('.')[0] for match in pdb_check(self.match_PDB_dir, base_only=True)]

        for match_name, match_df in df.groupby('match_name'):
            if match_name in match_PDBs:
                min_score = match_df.min()
                match_score_dict[match_name] = min_score['match_score']

        return match_score_dict

    def _import_ideal_binding_sites(self):
        """
        Import all ideal binding site PDBs whose constraints were used to find successful matches.
        Key: name of ideal binding motif without .pdb extension
        Value: prody object of ideal binding motif
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
        # print('Ligand 11A shell CB count: {}'.format(ligand_shell_eleven))

        # todo: UM_1_Y176W276Q170Q177_1_2BH1_TEP_0001_21-22-25-5_1.pdb doesn't have two chains???
        interface_CB_contact_percentage = 0
        motif_shell_CB = 0

        # Percentage of CB atoms in the protein-protein interface (based on an 8A° threshold) that are within 6A of any ligand atom
        chains_in_dimer = list(set(match_prody.getChids()) - set('X'))
        if len(chains_in_dimer) == 2:
            interface_cb = match_prody.select('(name CB and chain {}) within 8 of chain {} or\
             (name CB and chain {}) within 8 of chain {}'.format(chains_in_dimer[0], chains_in_dimer[1],
                                                                 chains_in_dimer[1], chains_in_dimer[0]))
            # print(len(interface_cb))

            ligand_shell_six = match_prody.select('name CB within 6 of resname {}'.format(self.ligand))
            interface_CB_contact_percentage = len(set(ligand_shell_six.getIndices()) & set(interface_cb.getIndices())) / len(interface_cb)
            # print('Interface CB contact percentage: {}'.format(interface_CB_contact_percentage))

        # neighbor bin of motif residues (i.e. number of CB atoms within 8A° of any motif residue CB atom)
        motif_resnums = [res[1] for res in motif_residue_IDs]

        # todo: update to accomodate arbitrary number of motif residues
        neighbor_bin_resnum_set = set()

        # for loop, get all CB within 8 of each motif residue
        for motif_resnum in motif_resnums:

            # set residue numbers
            motif_residue_shell_CB = match_prody.select('name CB within 8 of (resnum {})'.format(motif_resnum))
            motif_residue_shell_resnums = set(motif_residue_shell_CB.getResnums())
            neighbor_bin_resnum_set = neighbor_bin_resnum_set | motif_residue_shell_resnums

        # select name CB and resnum [list of residue numbers]
        motif_shell_CB = len(set(neighbor_bin_resnum_set))

        # print('Motif shell CB count: {}'.format(motif_shell_CB))
        return ligand_shell_eleven, interface_CB_contact_percentage, motif_shell_CB

    def calculate_rmsd_stats(self, match_prody, ideal_name, motif_residue_IDs, match_name):
        """
        RMSD things
        
        *IMPORTANT FOR DETERMINING IDEAL-MATCHED RESIDUE PAIRS*
        Ideal binding motifs and constraint files follow the same residue ordering
        Match residues in matched PDB names follow same order as constraint file
        
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

        # print('RMSD: {}'.format(prody.calcRMSD(transformed_ideal_prody.select('resname {}'.format(self.ligand)), match_atom_coords)))

        # Select residues from ideal and get coords

        hv = transformed_ideal_prody.getHierView()
        backbone_atom_name_list = ['N', 'CA', 'C', 'O', 'OXT']

        # ideal_residue_prody_list = [res.select('not hydrogen') for res in hv.iterResidues()]
        ideal_residue_prody_list = []
        backbone_contact_bool = []

        # make copy of ligand so atom indicies are reset
        ligand_copy = transformed_ideal_prody.select('resname {} and not hydrogen'.format(self.ligand)).copy()

        for motif_residue in hv.iterResidues():
            if motif_residue.getResname() != self.ligand:

                # Make copy of residue so atom indicies are reset
                ideal_motif_copy = motif_residue.select('not hydrogen').copy()

                # Get closest atom-atom contacts
                # row_index_low == motif_copy && column_index_low == ligand_copy
                contact_distance, row_index_low, column_index_low = minimum_contact_distance(ideal_motif_copy, ligand_copy, return_indices=True)

                # If backbone contact {C, CA, N, O} then only consider backbone atoms for RMSD calculation
                motif_contact_atom = ideal_motif_copy.select('index {}'.format(row_index_low))
                if motif_contact_atom.getNames()[0] in backbone_atom_name_list:
                    ideal_residue_prody_list.append(ideal_motif_copy.select('name {}'.format(' '.join(backbone_atom_name_list))))
                    backbone_contact_bool.append(True)
                else:
                    ideal_residue_prody_list.append(ideal_motif_copy.select('protein'))
                    backbone_contact_bool.append(False)

        # Select residues from match and get coords

        # motif_residue_prody_list = [match_prody.select('resnum {} and not hydrogen and protein'.format(res_tuple[1])) for res_tuple in motif_residue_IDs]
        motif_residue_prody_list = []

        # CANNOT USE MATCHED MOTIF RESIDUE TO GET CONTACTS, MATCHED MOTIF CAN BE ALL FUCKED UP
        for res_tuple, bb_contact in zip(motif_residue_IDs, backbone_contact_bool):
            # Make copy of residue so atom indicies are reset
            match_motif_copy = match_prody.select('resnum {} and not hydrogen and protein'.format(res_tuple[1])).copy()

            # # Get closest atom-atom contacts
            if bb_contact:
                motif_residue_prody_list.append(match_motif_copy.select('name {}'.format(' '.join(backbone_atom_name_list))))
            else:
                motif_residue_prody_list.append(match_motif_copy.select('protein'))

        # Calculate match score as defined by ME!
        # todo: UM_1_E252W248T285W6_1_3A4U_TEP_0001-11-17-2-22_1.pdb has trouble aligning...
        try:
            # residue_match_score = sum([prody.calcRMSD(ideal, match) for ideal, match in zip(ideal_residue_prody_list, motif_residue_prody_list)])
            residue_match_score = 0
            for ideal, match in zip(ideal_residue_prody_list, motif_residue_prody_list):
                print(ideal.getResnames()[0], len([a for a in ideal]), match.getResnames()[0], len([b for b in match]))
                residue_match_score += prody.calcRMSD(ideal, match)

        except Exception as e:
            print('\n\n{}\n\n'.format(e))
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

        # Consolidate gurobi solutions into a single dataframe for easy lookup
        # Pirated from motifs.Generate_Constraints
        gurobi_solutions = pd.DataFrame(columns=['Obj_score', 'Residue_indicies', 'Conformer'])

        for solution_set in os.listdir(args['<gurobi_solutions_dir>']):
            temp_solution_df = pd.read_csv(os.path.join(args['<gurobi_solutions_dir>'], solution_set), usecols=['Obj_score', 'Residue_indicies', 'Conformer'])
            gurobi_solutions = gurobi_solutions.append(temp_solution_df, ignore_index=True)

        # Calculate stats for each matched PDB
        def evaluate_match(matched_PDB):
            """
            Evaluates quality metrics for a given match
            :param matched_PDB: path to a single matched PDB
            """
            match_name = os.path.basename(os.path.normpath(matched_PDB))
            match_prody = prody.parsePDB(matched_PDB)

            # Parse matched_PDB to get ideal binding site name and residues
            match_pdb_name = os.path.basename(os.path.normpath(matched_PDB))
            pnc = re.split('_|-|\.', match_pdb_name)
            match_name_dash_split = re.split('-|\.', match_pdb_name)

            motif_index_list = match_name_dash_split[1].split('_')[1:-1]
            motif_index_string = '_'.join(motif_index_list)

            # Ideal Binding Site Name
            ideal_binding_site_name = '{}_{}-1_{}'.format(pnc[5], pnc[6], motif_index_string)

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
            # todo: update to accomodate arbitrary number of motif residues
            try:
                min_res_per_chain = min([motif_residue_chain_list.count(chain) for chain in (set(motif_residue_chain_list) - set('X'))])
                if min_res_per_chain == 4:
                    min_res_per_chain = 0
            except:
                min_res_per_chain = -1

            # Count CB that are within 2.4A of ligand
            clashing_CB_atoms = match_prody.select('name CB within 2.4 of chain X')
            ligand_CB_clashes = len(clashing_CB_atoms) if clashing_CB_atoms is not None else 0

            # Look up binding motif score in gurobi solutions
            current_conformer = '{}_{}'.format(pnc[5], pnc[6])
            index_list_string = '[1, {}]'.format(', '.join(motif_index_list))

            gurobi_score_row = gurobi_solutions.loc[(gurobi_solutions['Residue_indicies'] == index_list_string) & (gurobi_solutions['Conformer'] == current_conformer)]
            gurobi_score = gurobi_score_row['Obj_score'].iloc[0]

            # Aggragate results
            row_dict = {'match_name': match_name,
                        'ligand_shell_eleven': ligand_shell_eleven,
                        'interface_CB_contact_percentage': interface_CB_contact_percentage,
                        'motif_shell_CB': motif_shell_CB,
                        'residue_match_score': residue_match_score,
                        'ligand_match_score': ligand_match_score,
                        'min_res_per_chain': min_res_per_chain,
                        'gurobi_motif_score': gurobi_score,
                        'ligand_CB_clashes': ligand_CB_clashes
                        }

            return row_dict

        # Multiprocess match evaluation
        process = Pool()
        match_metrics_list_of_dicts = process.map(evaluate_match, [pdb for pdb in pdb_check(match_PDB_dir)])
        process.close()
        process.join()

        # Return passing matcher results
        df = pd.DataFrame(match_metrics_list_of_dicts)
        df.set_index(['match_name'], inplace=True)
        df.to_csv('Match_Filter_Results-{}.csv'.format(args['<ligand>']))

    # Let's say take top 5% of hits, for each metric, passing matcher results have to be in all top 5%
    if args['--csv']:
        df.set_index(['match_name'], inplace=True)
    percent_cutoff = 0.25
    cut_index = int(len(df) * percent_cutoff)

    # Percentage cutoff filters
    # Ascending = True or False
    ascending_dict = {'ligand_shell_eleven': False,
                      'interface_CB_contact_percentage': False,
                      'motif_shell_CB': False,
                      # 'residue_match_score': True,
                      # 'ligand_match_score': True,
                      # 'min_res_per_chain': True # Not necessary if iterating to build full binding sites
                      # 'gurobi_motif_score': True # I don't think we will want to filter based on motif scores...
                      }

    # Hard Cutoff filters
    # Accept anything up to value for each key
    hard_cutoff_dict = {'ligand_CB_clashes': 1}

    # Dump lists of passing matches for each metric into set_list
    set_list = []

    print('Evaluating matches based on the following criteria:')

    # Return top percentage of each metric as determined by {percent_cutoff}
    # Hard Cutoff filters as determined by {hard_cutoff_dict}
    for index, column in df.iteritems():
        if index in list(ascending_dict.keys()):
            print(index)
            set_list.append(set(column.sort_values(ascending=ascending_dict[index]).index.tolist()[:cut_index]))
        elif index in list(hard_cutoff_dict.keys()):
            print(index)
            set_list.append(set(column.loc[column <= hard_cutoff_dict[index]].index.tolist()))

    # Get intersection of matches that are in top percentage for each metric
    initial_set = set.intersection(*set_list)

    # todo: add option or setting for this filter based on iterative design steps... maybe we don't need this anymore?
    # I don't think we need this anymore...

    # # Get all matches with min_res_per_chain = 0
    # min_res_none_set = set(df.groupby(['min_res_per_chain']).get_group(0).index.tolist())
    #
    # # Remove matches from final pool if there are not motif residues shared by both partners
    # final_set = initial_set - min_res_none_set

    final_set = initial_set
    pprint.pprint(final_set)
    print(len(final_set))

    # Move matches that pass filters into a new directory
    filtered_matches_dir = 'Matches-Filtered'
    os.makedirs(filtered_matches_dir, exist_ok=True)
    for match in list(final_set):
        shutil.copy(os.path.join(match_PDB_dir, match), filtered_matches_dir)

    # Export csv with metrics of final set
    final_set_df = df.loc[list(final_set)]
    final_set_df.to_csv(os.path.join(filtered_matches_dir, 'Final_set_metrics.csv'))
