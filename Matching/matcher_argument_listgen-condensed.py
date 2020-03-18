#!/usr/bin/env python3

"""
Generate a list of argument lists for matcher cluster runs and export as json

{
    "target_compound": "REN",
    "posfile_dir": "/path/to/posfile_dir",
    "cst_dir": "path/to/cst_dir",
    "params_dir": "path/to/params_dir",
    "scaffold_dir": "path/to/scaffold_dir",
    "output_path": "path/to/output_path",
    "arg_list": [
                    {
                        "conformer": "REN_0001",
                        "cst_file": "REN_0001-1_2_3_4.cst",
                        "scaffold": "ASDF.pdb.gz",
                        "posfile": "ASDF.pdb.pos"
                    }

                ]
}

Usage:
    matcher_argument_listgen (monomer|dimer|recover) <target_compound> <working_dir_name> <scaffold_file_path>
    ( <cst_dir> | iteration <gurobi_constraints_csv> <gurobi_iteration_solutions>) [-m][-g][-j]

Arguments:
    monomer                 Match to monomer scaffold library
    dimer                   Match to dimer scaffold library
    recover                 Match to user defined scaffolds
    <target_compound>       Three letter code of target compound
    <working_dir_name>      Name of the directory under {bsff_path}/Compounds/{target_compound_code} on cluster
    <scaffold_file_path>    Path to file containing list of scaffold PDBs
    <cst_dir>               Path to directory containing the target constraint files
        
Options:
    -m --multiple_params    Multiple params files needed!
    -g --add_gridlig        Generate arguments with gridlig file specified. Gridlig files are not required by the matcher,
                            and quite frankly the way we are using them doesn't really help anything much...
    -j --james              Use the dimer scaffold library by james
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
bsff_path = os.path.join('/wynton', 'home', 'kortemme',  'james.lucas', 'BindingSitesFromFragments')
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

def add_arg_to_list(scaffold, constraint, monomer=False):
    """
    Add an arguement to the list

    dimer - ASDF.pdb.gz
    monomer - ASDF.pdb_0.pos

    :return:
    """
    # Get current scaffold
    current_scaffold = '{}.pdb'.format(scaffold.split('.')[0]) if monomer else '{0}'.format(scaffold)

    # Get posfile file for current scaffold
    posfile_name = scaffold if monomer else scaffold[:-3] + '.pos'

    # Get params file path
    conformer = constraint.split('-')[0]
    params_name = '{0}.params'.format(conformer)

    arg = [conformer, constraint, current_scaffold, posfile_name, params_name]

    if args['--add_gridlig']:
        gridlig_name = '{0}.gridlig'.format(scaffold[:-4]) if monomer else '{0}.gridlig'.format(current_scaffold)
        arg += gridlig_name

    argument_list.append(arg)

# Start json dict
json_dict = {}

json_dict["target_compound"] = target_compound_code
json_dict["scaffold_pdb_dir"] = scaffold_pdb_path
json_dict["cst_dir"] = constraint_file_path
json_dict["params_dir"] = os.path.join(target_compound_path, 'params')
json_dict["posfile_dir"] = scaffold_pdb_path if args['monomer'] else os.path.join(scaffold_path, 'posfiles')
json_dict["output_path"] = os.path.join(working_dir_path, 'matches')

if args['--add_gridlig']:
    json_dict['gridlig_path'] = os.path.join(scaffold_path, 'gridligs')

# Get scaffold file names (local)
scaffold_file = open(args['<scaffold_file_path>'], 'r')
scaffold_file_list = [file_name.strip() for file_name in scaffold_file]

if args['iteration']:
    print(args['<gurobi_constraints_csv>'])
    scaffold_csv = pd.read_csv(args['<gurobi_constraints_csv>'])
    gurobi_iteration_solutions_dir = args['<gurobi_iteration_solutions>']
    gurobi_iteration_solutions = [solution for solution in os.listdir(gurobi_iteration_solutions_dir) if solution.endswith('.csv')]

    relevant_scaffolds = [scaffold for scaffold in scaffold_file_list if any([pdb in scaffold for pdb in scaffold_csv['scaffold'].values])]

    for index, row in scaffold_csv.iterrows():
        current_solution_set = [a for a in gurobi_iteration_solutions if row['constraint'] in a]

        # todo: announce skipped consrtaints... this assumes all iteration solutions are being provided
        # assert len(current_solution_set) == 1

        if len(current_solution_set) < 1:
            print('Gurobi solutions for constraints {0} not provided. Moving on!'.format(row['constraint']))
            continue

        current_solution_csv = pd.read_csv(os.path.join(gurobi_iteration_solutions_dir, current_solution_set[0]), usecols=['Conformer', 'Obj_score', 'Residue_indicies'])
        current_scaffolds = [scaffold for scaffold in relevant_scaffolds if row['scaffold'] in scaffold]

        for scaffold in current_scaffolds:
            print('Writing constraints for scaffold {}'.format(scaffold))
            for index, row in current_solution_csv.iterrows():
                residue_index_list = literal_eval(row['Residue_indicies'])
                constraint_filename = "{}-{}.cst".format(row['Conformer'], '_'.join([str(a) for a in residue_index_list]))
                add_arg_to_list(scaffold, constraint_filename, monomer=args['monomer'])

else:
    print(args['<cst_dir>'])
    for scaffold in scaffold_file_list:
        for constraint in os.listdir(args['<cst_dir>']):
            add_arg_to_list(scaffold, constraint, monomer=args['monomer'])

json_dict['args'] = argument_list

##############################
#   Dump Arguments as JSON   #
##############################
print("Generated {} constraints".format(len(argument_list)))

with open('matcher_argument_list.json', 'w') as arg_list:
    json.dump(json_dict, arg_list)