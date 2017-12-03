#!/usr/bin/env python3

"""
Generate a list of argument lists for matcher cluster runs and export as json

Usage:
    matcher_argument_listgen (monomer|dimer|recover) <target_compound> <working_dir_name> <scaffold_file_path> ( <cst_dir> | iteration <gurobi_constraints_csv> <gurobi_iteration_solutions>) [-m][-g][-j]

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
import pandas as pd
from ast import literal_eval

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

################################
#   Construct Arguement List   #
################################
print(match_type_index)
print(scaffold_path)
print(scaffold_pdb_path)

# Start list of args
argument_list = []

def add_arg_to_list(scaffold, constraint):
    """
    Add an arguement to the list
    :return:
    """
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

# Get scaffold file names (local)
if args['<scaffold_file_path>']:
    scaffold_file = open(args['<scaffold_file_path>'], 'r')
    scaffold_file_list = [file_name.strip() for file_name in scaffold_file]

# Lazy
if args['iteration']:
    scaffold_csv = pd.read_csv(args['<gurobi_constraints_csv>'])
    gurobi_iteration_solutions_dir = args['<gurobi_iteration_solutions>']
    gurobi_iteration_solutions = [solution for solution in os.listdir(gurobi_iteration_solutions_dir) if solution.endswith('.csv')]

    relevant_scaffolds = [scaffold for scaffold in scaffold_file_list if any([pdb in scaffold for pdb in scaffold_csv['scaffold'].values])]

    for index, row in scaffold_csv.iterrows():
        current_solution_set = [a for a in gurobi_iteration_solutions if row['constraint'] in a]
        assert len(current_solution_set) == 1

        current_solution_csv = pd.read_csv(os.path.join(gurobi_iteration_solutions_dir, current_solution_set[0]), usecols=['Conformer', 'Obj_score', 'Residue_indicies'])
        current_scaffolds = [scaffold for scaffold in relevant_scaffolds if row['scaffold'] in scaffold]

        for scaffold in current_scaffolds:
            for index, row in current_solution_csv.iterrows():
                residue_index_list = literal_eval(row['Residue_indicies'])
                constraint_filename = "{}-{}.cst".format(row['Conformer'], '_'.join([str(a) for a in residue_index_list]))
                add_arg_to_list(scaffold, constraint_filename)

else:
    for scaffold in scaffold_file_list:
        for constraint in os.listdir(args['<cst_dir>']):
            add_arg_to_list(scaffold, constraint)

##############################
#   Dump Arguments as JSON   #
##############################
print("Generated {} constraints".format(len(argument_list)))

with open('matcher_argument_list.json', 'w') as arg_list:
    json.dump(argument_list, arg_list)