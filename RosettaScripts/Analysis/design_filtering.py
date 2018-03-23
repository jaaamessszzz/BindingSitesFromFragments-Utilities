#!/usr/bin/env python3

"""
Filter designs based on... stuff. Same idea as the match filtering script, populate a csv with quality metrics for
each design that I can use to rank the best designs.

NOTE: This script assumes that the Matcher output PDB is in the same directory as <design_json> (I know I'm going to
regret this later...)

Usage:
    design_filtering consolidate <design_PDB_dir> <target_designs> <design_json> [<score_term>...] [--fill_quota] [--unique] [--HBNet_911]
    design_filtering analyze <design_PDB_dir> <scores_csv> <design_json>

Arguments:
    consolidate
        Calculate and consolidate statistics for a design

    analyze
        Calculate more things

    <design_PDB_dir>
        Directory containing design PDBs and score.sc files

    <score_term>
        Score term(s) to filter on

    <target_designs>
        Number of designs return

Options:

    --fill_quota -f
        Keep filtering designs until <target_designs> are found

    --unique -u
        Return unique sequences only

    --HBNet_911 -h
        HBNet 911... parse pdb files to get correct HB unsat values

"""
import os
import sys
import re
import shutil
import pprint
import json

import docopt
import pandas as pd
import numpy as np
import prody

class FilterDesigns():
    """
    Modules for filtering designs
    """

    def __init__(self, design_dir, designs_to_return, update_hb_unsats, design_json_path):
        self.design_dir = design_dir
        self.design_json = json.load(open(design_json_path, 'r'))
        self.match_input_path = os.path.join(os.path.dirname(design_json_path), self.design_json['pdbid'])
        self.update_hb_unsats = update_hb_unsats
        self.df = self.consolidate_scores()
        self.designs_to_return = designs_to_return
        self.processed_files = set()
        self.skipped_files = set()
        self.returned_sequences = set()
        self.unique_sequences = set()

    def copy_designs(self, design_list, directory_name='Fresh_Designs', unique_only=False):
        """
        Copies the specified PDBs to new directory

        :param directory_name: Destination directory of copied designs
        :param design_list: list of PDBs (without .pdb extension)
        :return:
        """

        print('Attempting to copy {0} designs...'.format(len(design_list)))

        # Make output directory in cwd
        dst_dir = os.path.join(os.getcwd(), directory_name)
        os.makedirs(directory_name, exist_ok=True)

        # Keep track of processed files
        current_designs_processed = set()
        current_designs_skipped = set()

        # Walk through design_dir to find all designs in design_list
        for dir_path, dir_name, file_name in os.walk(self.design_dir):
            if directory_name not in dir_path:
                files_to_transfer = set(file_name) & set(design_list)

                if unique_only:
                    for file in files_to_transfer:
                        # Unique sequence check
                        hv = prody.parsePDB(os.path.join(dir_path, file)).getHierView()
                        seq_tuple = tuple([chain.getSequence(allres=True) for chain in hv])

                        # If PDB sequence is unique, copy to dst
                        if seq_tuple not in self.unique_sequences:
                            self.unique_sequences.add(seq_tuple)
                            shutil.copy2(os.path.join(dir_path, file), dst_dir)
                            current_designs_processed.add(file)
                        else:
                            current_designs_skipped.add(file)
                            print('{0} skipped: PDB with identical sequence already identified'.format(file))

                else:
                    for file in files_to_transfer:
                        shutil.copy2(os.path.join(dir_path, file), dst_dir)
                        current_designs_processed.add(file)

        # Report what files are missing
        missing_files = set(design_list) - (current_designs_processed | current_designs_skipped)
        if len(missing_files) == 0 and len(current_designs_skipped) == 0:
            print('\nAll files were transferred successfully!')

        if len(missing_files) != 0:
            print('\nThe following {0} files were not found in {1}:'.format(len(missing_files), self.design_dir))
            for file in missing_files:
                print(file)

        if len(current_designs_skipped) !=0:
            print('\nThe following {0} files were skipped:'.format(len(current_designs_skipped)))
            for file in current_designs_skipped:
                print(file)

        for design in current_designs_processed:
            self.processed_files.add(design)
        for design in current_designs_skipped:
            self.skipped_files.add(design)

    def get_best_designs(self, filter_terms):
        """
        Get designs that are in the top XX% for all desired metrics
        :param df: dataframe containing scoring information for designs
        :return:
        """

        # todo: update to be more generalizable... only takes 50 best scoring designs at the moment...
        # default = value at which to filter if design is above/below
        # ascending = True -> Low to high, False -> High to low
        recognized_filter_terms = {'total_score': {'ascending': True, 'default': 0},
                                   'Lig_ShapeC': {'ascending': False, 'default': 0.55},
                                   'buried_unsat': {'ascending': True, 'default': 30},
                                   'buried_unsat_bb': {'ascending': True, 'default': 10},
                                   'buried_unsat_sc': {'ascending': True, 'default': 10},
                                   'interfE': {'ascending': True, 'default': 0},
                                   'packstat': {'ascending': False, 'default': 0.5},
                                   'rosettaholes': {'ascending': True, 'default': 1},
                                   'ligand_fa_rep': {'ascending': True, 'default': 10},
                                   'ligand_total': {'ascending': True, 'default': 0}
                                   }
        set_list = []

        # Default to total_score if filter_terms == []
        if filter_terms == []:
            filter_terms = ['total_score']

        # Iterate through dataframe and add list of passing designs for each criterion to set
        for index, column in self.df.iteritems():
            current_set = set()

            if index in filter_terms:
                # sorted_column = column.sort_values(ascending=recognized_filter_terms[index]['ascending']).index.tolist()
                sorted_column = column.sort_values(ascending=recognized_filter_terms[index]['ascending'])

                if recognized_filter_terms[index]['ascending'] is True:
                    filtered_sorted_column = sorted_column[sorted_column < recognized_filter_terms[index]['default']].index.tolist()
                else:
                    filtered_sorted_column = sorted_column[sorted_column > recognized_filter_terms[index]['default']].index.tolist()

                print(index)
                print(len(filtered_sorted_column))

                # todo: reintegrate this somehow... <target_designs> is not respected at the moment...
                # added_to_set = []
                #
                # for value in filtered_sorted_column:
                #     if '{0}.pdb'.format(value) not in self.processed_files:
                #         current_set.add('{0}.pdb'.format(value))
                #         self.processed_files.add('{0}.pdb'.format(value))
                #         added_to_set.append('{0}.pdb'.format(value))
                #
                # if len(added_to_set) == self.designs_to_return:
                #     set_list.append(current_set)
                #     break

                set_list.append(set(filtered_sorted_column))

        # Take intersection of all filter metrics so we only end up with designs that pass all specified filters
        initial_set = set.intersection(*set_list)

        # --- Pirated from match_filtering.py --- #
        # Only return designs that score in the top 25th percentile of all score terms specified

        df_subset = self.df.loc[list(initial_set)]
        cut_index = int(len(df_subset) * 0.25)
        percentage_set_list = []

        for index, column in df_subset.iteritems():
            if index in filter_terms:
                percentage_set_list.append(set(column.sort_values(ascending=recognized_filter_terms[index]['ascending']).index.tolist()[:cut_index]))

        # Get intersection of matches that are in top percentage for each metric
        final_set = set.intersection(*percentage_set_list)

        return ['{0}.pdb'.format(pdb) for pdb in final_set]

    def convert_score_to_dict(self, score_file):
        """
        Converts score files into a dict

        Returns list of dicts (thinking ahead...)
        :return:
        """
        converted_dict_list = []
        score_terms = []
        score_values_list = []

        # Populate score_terms and score_values_list with values from score file
        with open(score_file, 'r') as file:
            for line in file:
                split_line = line.split()

                # Score terms
                if split_line[0].strip() == 'SCORE:' and split_line[1].strip() == 'total_score':
                    score_terms = [term.strip() for term in split_line[1:]]

                # Score values
                elif split_line[0].strip() == 'SCORE:':
                    score_values_line = [float(term.strip()) for term in split_line[1:-1]]
                    score_values_line.append(split_line[-1])
                    score_values_list.append(score_values_line)

        # Convert score_values_list into a list of dicts with keys score_terms
        for score_list in score_values_list:
            converted_dict_list.append({score_term: score_value for score_term, score_value in zip(score_terms, score_list)})

        return converted_dict_list

    def consolidate_scores(self):
        """
        Consolidate score files into a dataframe that can be exported as a csv
        :return:
        """
        # Store all scores in dict of dicts
        temp_dict = {}

        # todo: convert to os.walk() to evaluate nested directories, all jobs currently submitted in the same directory
        for file in os.listdir(self.design_dir):
            if file.endswith('.sc'):

                # Convert scorefile to dict
                converted_scorefile_dict = self.convert_score_to_dict(os.path.join(self.design_dir, file))

                # Add to scores to dict
                for score_dict in converted_scorefile_dict:
                    temp_dict[score_dict['description']] = score_dict

        # Convert dict into dataframe
        df = pd.DataFrame.from_dict(temp_dict, orient='index')

        # Input PDB name from output
        template_name = df['description'][0]
        temp_1 = template_name.split('.')[0]
        PDB_basename = re.split('-', temp_1, maxsplit=1)[1][:-5]

        # Set index
        df.set_index('description', inplace=True)

        # --- Process HBUnsats and other metrics written to PDB --- #
        # I've decided this is too much of a pain in the ass to implement... organizationally...

        # # Get HB Unsat resnums from WT scaffold PDB
        # with open(self.match_input_path) as derp:
        #     for line in derp:
        #         if line.startswith(('ATOM', 'HETATM', 'CONECT', 'TER')):
        #             continue
        #
        #         if line.startswith('BuriedUnsatHbonds buried_unsat:'):
        #             existing_unsats = set()
        #             continue
        #
        #         if line.startswith('      Unsatisfied H'):
        #             existing_unsats.add(int(re.split('[ :]', line)[6]))
        #             continue
        #
        #         if line.startswith('# All scores below are weighted scores, not raw scores.'):
        #             break

        # Get residues within 5A of design/motif residues that can potentially have unsats from designs
        input_PDB_prody = prody.parsePDB(self.match_input_path)

        designable_and_motif_positions = [a[1] for a in self.design_json['match_residues']] + self.design_json['design_residue_list']
        potential_unsat_set = set(designable_and_motif_positions)

        for residue in designable_and_motif_positions:
            ca_resnum_list = [res.getResnum() for res in input_PDB_prody.select('protein within 5 of resnum {}'.format(residue))]
            potential_unsat_set.update(ca_resnum_list)

        # Please forgive this nested mess...
        for file in os.listdir(self.design_dir):
            file_path = os.path.join(self.design_dir, file)

            # Go though all pdb files...
            if file.endswith('.pdb'):

                final_bb_unsat_value = final_sc_unsat_value = ligand_total = ligand_fa_rep = 0

                # Save heavy unsat count
                with open(file_path, 'r') as sf:
                    ligand_id = None

                    for line in sf:
                        if line.startswith(('ATOM', 'HETATM', 'CONECT', 'TER')):
                            continue

                        # --- Process HBNet Blocks --- #

                        # HBUnsats gets written several times... the "real" value is the last one written to file
                        if line.startswith('BuriedUnsatHbonds buried_unsat:'):
                            final_bb_unsat_value = final_sc_unsat_value = 0
                            relevant_unsat_set = set()
                            continue

                        if line.startswith('      Unsatisfied H'):
                            split_line = [a for a in re.split('[ :]', line.strip()) if a]

                            # Only count unsats if they are in potential_unsat_set
                            if int(split_line[6]) in potential_unsat_set:
                                if 'H' in split_line[-1] and split_line[1] == 'Hpol':
                                    continue

                                elif split_line[-1] in ['N', 'O']:
                                    final_bb_unsat_value += 1

                                else:
                                    final_sc_unsat_value += 1

                                # Keep track of resnums for counted unsats
                                relevant_unsat_set.add(split_line[6])

                        if line.startswith('{}_'.format(ligand_id)) and ligand_id is not None:
                            ligand_fa_rep = float(line.split()[2])
                            ligand_total = float(line.split()[20])

                        # Get ligand code
                        if line.startswith('HETNAM'):
                            ligand_id = line.split()[1]
                            continue

                # --- Fix values in df --- #
                df.at[file.split('.')[0], 'buried_unsat'] = final_bb_unsat_value + final_sc_unsat_value
                df.at[file.split('.')[0], 'buried_unsat_bb'] = final_bb_unsat_value
                df.at[file.split('.')[0], 'buried_unsat_sc'] = final_sc_unsat_value
                df.at[file.split('.')[0], 'relevant_unsats'] = '+'.join(list(relevant_unsat_set))
                df.at[file.split('.')[0], 'ligand_fa_rep'] = ligand_fa_rep
                df.at[file.split('.')[0], 'ligand_total'] = ligand_total

        df.to_csv('{}-Consolidated_Scores.csv'.format(os.path.basename(self.design_dir)))

        return df


class AnalyzeDesigns():

    def __init__(self, input_design_dir, csv_path, design_json_path):
        self.design_dir = input_design_dir
        self.df = pd.read_csv(csv_path)
        self.design_json = json.load(open(design_json_path, 'r'))

    def sequence_logo(self):
        """
        Generate sequence logo for designed positions
        :return:
        """

        # Alternatively, list of lists for each design
        design_lol = []

        for file in os.listdir(self.design_dir):
            if file.endswith('.pdb'):

                temp_list = []

                hv = prody.parsePDB(os.path.join(self.design_dir, file)).getHierView()
                seq = ''.join([chain.getSequence(allres=True) for chain in hv])

                for design_position in self.design_json['design_residue_list']:
                    temp_list.append(seq[design_position - 1])

                design_lol.append(''.join(temp_list))

        from Bio.Seq import Seq
        from Bio import motifs
        from Bio.Alphabet import IUPAC

        motif = motifs.create([Seq(designed_positions, IUPAC.protein) for designed_positions in design_lol])
        motif.weblogo("{0}-{1}-Weblogo.pdf".format(self.design_json['pdbid'].split()[0], os.path.basename(os.path.normpath(self.design_dir))), format='PDF')


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    input_design_dir = args['<design_PDB_dir>']
    filter_score_terms = args['<score_term>']
    fill_quota = args['--fill_quota']
    return_unique_only = args['--unique']
    update_hb_unsats = args['--HBNet_911']
    design_json_path = args['<design_json>']

    if args['consolidate']:
        target_design_number = int(args['<target_designs>'])

        filter_designs = FilterDesigns(input_design_dir, target_design_number, update_hb_unsats, design_json_path)

        if fill_quota:
            while len(filter_designs.unique_sequences) < target_design_number and len(filter_designs.processed_files) < len(filter_designs.df):
                design_list = filter_designs.get_best_designs(filter_score_terms)
                filter_designs.copy_designs(design_list, unique_only=return_unique_only)
                print('\n\n{0} unique designs have been found so far...\n\n'.format(len(filter_designs.unique_sequences)))
            print('\n\n{0} designs have been evaluated to return {1} unique designs!\n\n'.format(len(filter_designs.processed_files), len(filter_designs.unique_sequences)))

        else:
            design_list = filter_designs.get_best_designs(filter_score_terms)
            filter_designs.copy_designs(design_list, unique_only=return_unique_only)

            # Return csv with metrics for design_list subset
            design_list_df = filter_designs.df.loc[[pdb[:-4] for pdb in design_list]]
            design_list_df.to_csv('{}-Consolidated_Scores-Design_List_Subset.csv'.format(os.path.basename(input_design_dir)))

    if args['analyze']:
        csv_path = args['<scores_csv>']
        analyze_designs = AnalyzeDesigns(input_design_dir, csv_path, design_json_path)
        analyze_designs.sequence_logo()