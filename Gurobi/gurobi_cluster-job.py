import sys
import os
import sqlite3
import subprocess
import pandas as pd

#################
# Set variables #
#################

print 'Starting Gurobi...'

# Add environment variables
os.environ['PATH'] += os.pathsep + '/netapp/home/james.lucas/gurobi751/linux64/bin'
os.environ['LD_LIBRARY_PATH'] += os.pathsep + '/netapp/home/james.lucas/gurobi751/linux64/lib'

os.environ['GUROBI_HOME'] = os.pathsep + '/netapp/home/james.lucas/gurobi751/linux64'
os.environ['GRB_LICENSE_FILE'] = os.pathsep + '/netapp/home/james.lucas/gurobi751/linux64/gurobi.lic'

# DEBUGGING
get_env = subprocess.Popen(['env'])
get_env.wait()

from gurobipy import *

struct_id = int(sys.argv[1])

#####################################################
# Select pairwise scores for structid aka conformer #
#####################################################

# Generate gurobi input table/csv
connection = sqlite3.connect('two_body_terms.db')
cursor = connection.cursor()

residue_table = pd.read_sql_query("SELECT * from residues", connection)
score_table = pd.read_sql_query(
    """
    SELECT struct_id, resNum1, resNum2, 
    CASE
    WHEN sum(score_value) > 0 THEN 1
    ELSE sum(score_value)
    END 
    as score_total from relevant_2b_scores where struct_id={0} group by struct_id, resNum1, resNum2;
    """.format(struct_id), connection)

#######################
# Solve for solutions #
#######################

# Set up model
residue_interactions = Model("residue_interactions")

# Add MIP binary variables
MIP_var_list = []
for index, row in residue_table.iterrows():
    MIP_var_list.append(residue_interactions.addVar(vtype=GRB.BINARY, name=str(row['resNum'])))

# Set up dict with pairwise scores
score_dict = {}
for index, row in score_table.iterrows():
    score_dict[(row['resNum1'], row['resNum2'])] = row['score_total']

# Set objective function
# Get all possible reisude pairs
residue_pairs = itertools.combinations(MIP_var_list, 2)
residue_interactions.setObjective(quicksum((MIP_var_list[int(key[0] - 1)] * MIP_var_list[int(key[1] - 1)] * value) for key, value in score_dict.items()), GRB.MINIMIZE)

# Add constraints
residue_interactions.addConstr(quicksum(MIP_var_list) <= 10) # Number of residues in a binding motif (includes ligand)
residue_interactions.addConstr(MIP_var_list[0] == 1)
for index, row in score_table.iterrows():
    if row['score_total'] > 0:
        residue_interactions.addConstr(MIP_var_list[int(row['resNum1'] - 1)] + MIP_var_list[int(row['resNum2'] - 1)] <= 1)

# Set Parameters
residue_interactions.Params.PoolSolutions = 10000
residue_interactions.Params.PoolGap = 0.2
residue_interactions.Params.PoolSearchMode = 2
residue_interactions.Params.Threads = 28

# Optimize
residue_interactions.optimize()

######################
# Retrieve solutions #
######################

# Get residue-index mapping for current conformer
index_mapping = pd.read_sql_query("SELECT * from residue_index_mapping WHERE struct_id = {0}".format(struct_id), connection, index_col='residue_index')

results_list = []

for i in range(residue_interactions.SolCount):
    residue_interactions.setParam(GRB.Param.SolutionNumber, i)

    res_index_tuple = [index + 1 for index, value in enumerate(residue_interactions.getVars()) if int(value.Xn) == 1]
    source_pdb_list = [index_mapping.loc[idx, 'source_pdb'] for idx in res_index_tuple if idx != 1]
    print('Residue Indicies: {}'.format(res_index_tuple))
    print('Obj: {}'.format(residue_interactions.objVal))
    print('Non-Ideal Obj: {}'.format(residue_interactions.PoolObjVal))

    results_list.append({'Residue_indicies': res_index_tuple,
                         'Obj_score': residue_interactions.PoolObjVal})

df = pd.DataFrame(results_list)
df.to_csv('Gurobi_results-{0}.csv'.format(struct_id))