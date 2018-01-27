#!/usr/bin/env python3

"""
Filter designs based on... stuff. Same idea as the match filtering script, populate a csv with quality metrics for
each design that I can use to rank the best designs.

Usage:
    design_filtering <design_PDB_dir>  <target_designs> [<score_term>...] [--fill_quota]

Arguments:
    <design_PDB_dir>
        Directory containing design PDBs and score.sc files

    <score_term>
        Score term(s) to filter on

    <target_designs>
        Number of designs return

Options:

    --fill_quota -f
        Keep filtering designs until

"""
import os
import sys
import re
import shutil
import pprint

import docopt
import pandas as pd
import numpy as np
import prody

class FilterDesigns():
    """
    Modules for filtering designs
    """

    def __init__(self, design_dir, designs_to_return):
        self.design_dir = design_dir
        self.df = self.consolidate_scores()
        self.designs_to_return = designs_to_return
        self.processed_files = set()
        self.skipped_files = set()
        self.unique_sequences = set()

    def copy_designs(self, design_list, directory_name='Fresh_Designs'):
        """
        Copies the specified PDBs to new directory

        :param design_list: list of PDBs (without .pdb extension)
        :return:
        """

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
                    if len(self.unique_sequences) == self.designs_to_return:
                        break

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
        recognized_filter_terms = {'total_score': False}
        set_list = []

        # Default to total_score if filter_terms == []
        if filter_terms == []:
            filter_terms = ['total_score']

        # Iterate through dataframe and add list of passing designs for each criterion to set
        for index, column in self.df.iteritems():
            current_set = set()

            if index in filter_terms:
                sorted_column = column.sort_values(ascending=recognized_filter_terms[index]).index.tolist()
                added_to_set = []

                for value in sorted_column:
                    if '{0}.pdb'.format(value) not in self.processed_files:
                        current_set.add('{0}.pdb'.format(value))
                        self.processed_files.add('{0}.pdb'.format(value))
                        added_to_set.append('{0}.pdb'.format(value))

                    if len(added_to_set) == self.designs_to_return:
                        set_list.append(current_set)
                        break

        # Take intersection of all filter metrics so we only end up with designs that pass all specified filters
        final_set = set.intersection(*set_list)

        return list(final_set)

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
                converted_scorefile_dict = self.convert_score_to_dict(file)

                # Add to scores to dict
                for score_dict in converted_scorefile_dict:
                    temp_dict[score_dict['description']] = score_dict

        # Convert dict into dataframe
        df = pd.DataFrame.from_dict(temp_dict, orient='index')

        return df


class AnalyzeDesigns():

    def __init__(self):
        self.herp = True

    def sequence_logo(self):
        """
        Generate sequence logo for designed positions
        :return:
        """
        return 0


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    input_design_dir = args['<design_PDB_dir>']
    filter_score_terms = args['<score_term>']
    target_design_number = int(args['<target_designs>'])
    fill_quota = args['--fill_quota']

    filter_designs = FilterDesigns(input_design_dir, target_design_number)

    if fill_quota:
        while len(filter_designs.unique_sequences) < target_design_number and len(filter_designs.processed_files) < len(filter_designs.df):
            design_list = filter_designs.get_best_designs(filter_score_terms)
            filter_designs.copy_designs(design_list)
            print('\n\n{0} unique designs have been found so far...\n\n'.format(len(filter_designs.unique_sequences)))
        print('\n\n{0} designs have been evaluated to return {1} unique designs!\n\n'.format(len(filter_designs.processed_files), len(filter_designs.unique_sequences)))

    else:
        design_list = filter_designs.get_best_designs(filter_score_terms, target_design_number)
        filter_designs.copy_designs(design_list)