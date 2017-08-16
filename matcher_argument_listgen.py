#!/usr/bin/env python3

"""
Generate a list of argument lists for matcher cluster runs and export as json

Usage:
    matcher_argument_listgen <target_compound> <working_dir_name> <scaffold_file_path> <cst_dir>

Arguments:     
    <target_compound>
        Three letter code of target compound
    
    <working_dir_name> 
        Name of the directory under {bsff_path}/Compounds/{target_compound_code} on cluster
        
    <scaffold_file_path>
        Path to file containing list of scaffold PDBs
    
    <cst_dir>
        Path to directory containing the target constraint files
        
Options:
        
"""
import os
import sys
import json
import docopt

args = docopt.docopt(__doc__)

###########################
# Variable and Path Setup #
###########################

# Get Target Compound Directory as argv
target_compound_code = args['<target_compound>']

# Define relevant paths (Cluster paths!!!)
bsff_path = os.path.join('/netapp', 'home', 'james.lucas', 'BindingSitesFromFragments')
scaffold_path = os.path.join(bsff_path, 'Dimer_Scaffolds')
target_compound_path = os.path.join(bsff_path, 'Compounds', target_compound_code)
working_dir_path = os.path.join(target_compound_path, args['<working_dir_name>'])
constraint_file_path = os.path.join(working_dir_path, 'cst_files')

# Get number of constraint files for target compound (local)
number_of_cst_files = len(os.listdir(args['<cst_dir>']))

# Get scaffold file names (local)
scaffold_file = open(args['<scaffold_file_path>'], 'r')
scaffold_file_list = [file_name.strip() for file_name in scaffold_file]

################################
#   Construct Arguement List   #
################################

# Start list of args
argument_list = []

for scaffold in scaffold_file_list:
    for constraint in os.listdir(args['<cst_dir>']):
        # Get current scaffold
        current_scaffold_path = os.path.join(scaffold_path, 'cleaned_heterodimers_all_biological_units', scaffold)

        # Get posfile and gridlig files for current scaffold
        posfile_name = scaffold[:-3] + '.pos'
        posfile_path = os.path.join(scaffold_path, 'posfiles', posfile_name)

        gridlig_name = scaffold[:-3] + '.gridlig'
        gridlig_path = os.path.join(scaffold_path, 'gridligs', gridlig_name)

        # Get params file path
        params_name = target_compound_code + '.params'
        params_path = os.path.join(target_compound_path, params_name)

        # Output path
        output_path = os.path.join(working_dir_path, 'matches')

        arg = [current_scaffold_path,
               target_compound_code,
               gridlig_path,
               posfile_path,
               os.path.join(constraint_file_path, constraint),
               params_path,
               output_path
               ]

        argument_list.append(arg)

with open('matcher_argument_list.json', 'w') as arg_list:
    json.dump(argument_list, arg_list)