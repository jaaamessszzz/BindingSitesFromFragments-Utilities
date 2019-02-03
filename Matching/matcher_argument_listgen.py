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

if args['dimer']:
    scaffold_path = os.path.join(bsff_path, 'Dimer_Scaffolds_by_James') if args['--james'] \
        else os.path.join(bsff_path, 'Dimer_Scaffolds')
    scaffold_pdb_path = os.path.join(scaffold_path, 'PDBs') if args['--james'] \
        else os.path.join(scaffold_path, 'cleaned_heterodimers_all_biological_units')

elif args['monomer']:
    scaffold_path = scaffold_pdb_path = os.path.join(bsff_path, 'Monomer_Scaffolds')

elif args['recover']:
    scaffold_path = os.path.join(working_dir_path, 'Recovery_Scaffolds')
    scaffold_pdb_path = os.path.join(scaffold_path, 'PDBs')
else:
    sys.exit('wut.')

################################
#   Construct Arguement List   #
################################
print(scaffold_path)
print(scaffold_pdb_path)

# Start list of args
argument_list = []

def add_arg_to_list_dimers(scaffold, constraint):
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

def add_arg_to_list_monomers(posfile, constraint):
    """
    Add a monomer scaffold to the list
    :param scaffold:
    :param constraint:
    :return:
    """
    scaffold_base = posfile.split('.')[0]

    # Get current scaffold
    scaffold = '{}.pdb'.format(scaffold_base)
    current_scaffold_path = os.path.join(scaffold_pdb_path, scaffold)

    # Get posfile and gridlig files for current scaffold
    posfile_path = os.path.join(scaffold_path, posfile)
    gridlig_path = os.path.join(scaffold_path, '{}.gridlig'.format(scaffold_base))

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
scaffold_file = open(args['<scaffold_file_path>'], 'r')
scaffold_file_list = [file_name.strip() for file_name in scaffold_file]

# Lazy
if args['dimer'] or args['recover']:

    if args['iteration']:
        print(args['<gurobi_constraints_csv>'])
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
                print('Writing constraints for scaffold {}'.format(scaffold))
                for index, row in current_solution_csv.iterrows():
                    residue_index_list = literal_eval(row['Residue_indicies'])
                    constraint_filename = "{}-{}.cst".format(row['Conformer'], '_'.join([str(a) for a in residue_index_list]))
                    add_arg_to_list_dimers(scaffold, constraint_filename)

    else:
        for scaffold in scaffold_file_list:
            for constraint in os.listdir(args['<cst_dir>']):
                add_arg_to_list_dimers(scaffold, constraint)
# Super lazy
elif args['monomer'] :
    # scaffold_file_list for monomers is a list of posfiles, pdb scaffolds in clean.list have suffix in *_11.pdb.clean
    # This is a little messed up since the EnzDes scaffold library has multiple posfiles for each scaffold...
    if args['iteration']:
        print(args['<gurobi_constraints_csv>'])
        scaffold_csv = pd.read_csv(args['<gurobi_constraints_csv>'])
        gurobi_iteration_solutions_dir = args['<gurobi_iteration_solutions>']
        gurobi_iteration_solutions = [solution for solution in os.listdir(gurobi_iteration_solutions_dir) if solution.endswith('.csv')]

        relevant_scaffolds = [scaffold for scaffold in scaffold_file_list if any([pdb in scaffold for pdb in scaffold_csv['scaffold'].values])]
        print(scaffold_file_list)
        print(relevant_scaffolds)

        for index, row in scaffold_csv.iterrows():
            current_solution_set = [a for a in gurobi_iteration_solutions if row['constraint'] in a]
            assert len(current_solution_set) == 1

            current_solution_csv = pd.read_csv(os.path.join(gurobi_iteration_solutions_dir, current_solution_set[0]), usecols=['Conformer', 'Obj_score', 'Residue_indicies'])
            current_scaffolds = [scaffold for scaffold in relevant_scaffolds if row['scaffold'] in scaffold]

            for scaffold in current_scaffolds:
                print('Writing constraints for scaffold {}'.format(scaffold))
                for index, row in current_solution_csv.iterrows():
                    residue_index_list = literal_eval(row['Residue_indicies'])
                    constraint_filename = "{}-{}.cst".format(row['Conformer'], '_'.join([str(a) for a in residue_index_list]))

                    # Required due to janky naming system and how I parse results after match filtering...
                    # scaffold_actual = scaffold if len(scaffold) == 4 else '{0}_11'.format(scaffold)

                    add_arg_to_list_monomers(scaffold, constraint_filename)

    else:
        for scaffold in scaffold_file_list:
            for constraint in os.listdir(args['<cst_dir>']):
                add_arg_to_list_monomers(scaffold, constraint)

##############################
#   Dump Arguments as JSON   #
##############################
print("Generated {} constraints".format(len(argument_list)))

with open('matcher_argument_list.json', 'w') as arg_list:
    json.dump(argument_list, arg_list)