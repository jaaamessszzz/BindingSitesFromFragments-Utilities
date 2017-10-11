#!/usr/bin/env python3

"""
Generate a list of argument lists for matcher cluster runs and export as json

Usage:
    matcher_argument_listgen (monomer|dimer|recover) <target_compound> <working_dir_name> <scaffold_file_path> <cst_dir> [-m][-g][-j]

Arguments:
    monomer
        Match to monomer scaffold library
        
    dimer
        Match to dimer scaffold library
         
     recover
        Match to user defined scaffolds
        
    <target_compound>
        Three letter code of target compound
    
    <working_dir_name> 
        Name of the directory under {bsff_path}/Compounds/{target_compound_code} on cluster
        
    <scaffold_file_path>
        Path to file containing list of scaffold PDBs
    
    <cst_dir>
        Path to directory containing the target constraint files
        
Options:
    -m --multiple_params
        Multiple params files needed!
        
    -g --add_gridlig
        Generate arguments with gridlig file specified. Gridlig files are not required by the matcher, and quite frankly 
        the way we are using them doesn't really help anything much...
        
    -j --james
        Use the dimer scaffold library by james
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
target_compound_path = os.path.join(bsff_path, 'Compounds', target_compound_code)
working_dir_path = os.path.join(target_compound_path, 'Matching', args['<working_dir_name>'])
constraint_file_path = os.path.join(working_dir_path, 'cst_files')

match_type_list = [args['monomer'], args['dimer'], args['recover']]
match_type_index = [i for i, x in enumerate(match_type_list) if x][0]

if match_type_index in [0, 1]:
    scaffold_dirname_list = ['Monomer_Scaffolds', 'Dimer_Scaffolds']
    scaffold_pdb_dirname_list = ['cleaned_heterodimers_all_biological_units', 'cleaned_heterodimers_all_biological_units']
    
    if args['--james'] and match_type_index == 1:
        scaffold_path = os.path.join(bsff_path, 'Dimer_Scaffolds_by_James')
        scaffold_pdb_path = os.path.join(scaffold_path, 'PDBs')
    else:
        scaffold_path = os.path.join(bsff_path, scaffold_dirname_list[match_type_index])
        scaffold_pdb_path = os.path.join(scaffold_path, scaffold_pdb_dirname_list[match_type_index])
        
else:
    scaffold_path = os.path.join(working_dir_path, 'Recovery_Scaffolds')
    scaffold_pdb_path = os.path.join(scaffold_path, 'PDBs')

# Get number of constraint files for target compound (local)
number_of_cst_files = len(os.listdir(args['<cst_dir>']))

# Get scaffold file names (local)
scaffold_file = open(args['<scaffold_file_path>'], 'r')
scaffold_file_list = [file_name.strip() for file_name in scaffold_file]

################################
#   Construct Arguement List   #
################################
print(match_type_index)
print(scaffold_path)
print(scaffold_pdb_path)

# Start list of args
argument_list = []

for scaffold in scaffold_file_list:
    for constraint in os.listdir(args['<cst_dir>']):
        # Get current scaffold
        current_scaffold_path = os.path.join(scaffold_pdb_path, scaffold)

        # Get posfile and gridlig files for current scaffold
        # todo: add check for compressed files
        # [:-3] lazily removes .gz suffix
        posfile_name = scaffold[:-3] + '.pos'
        posfile_path = os.path.join(scaffold_path, 'posfiles', posfile_name)

        gridlig_name = scaffold[:-3] + '.gridlig'
        gridlig_path = os.path.join(scaffold_path, 'gridligs', gridlig_name)

        # Get params file path
        if args['--multiple_params']:
            params_name = os.path.join('params', constraint.split('-')[0] + '.params')
        else:
            params_name = target_compound_code + '.params'

        params_path = os.path.join(target_compound_path, params_name)

        # Output path
        output_path = os.path.join(working_dir_path, 'matches')

        arg = [current_scaffold_path,
               target_compound_code,
               posfile_path,
               os.path.join(constraint_file_path, constraint),
               params_path,
               output_path
               ]

        if args['--add_gridlig']:
            arg.append(gridlig_path)

        argument_list.append(arg)

with open('matcher_argument_list.json', 'w') as arg_list:
    json.dump(argument_list, arg_list)