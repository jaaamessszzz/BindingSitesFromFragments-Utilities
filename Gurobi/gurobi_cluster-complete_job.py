import socket
from datetime import *
import sys
import os
import sqlite3
import subprocess
import pandas as pd
import shutil
import json
import time
import re

print("Python version:", sys.version)
print("Hostname:", socket.gethostname())

######################################################
# Timing things from Kyle, borrowed way back when... #
######################################################

def roundTime(dt=None, roundTo=1):
    """
    Round a datetime object to any time period (in seconds)
    dt : datetime.datetime object, default now.
    roundTo : Closest number of seconds to round to, default 1 second.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    http://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python/10854034#10854034
    """
    if dt == None : dt = datetime.now()
    seconds = total_seconds(dt - dt.min)
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + timedelta(0,rounding-seconds,-dt.microsecond)

def total_seconds(td):
    '''
    Included in python 2.7 but here for backwards-compatibility for old Python versions (like on QB3 cluster)
    '''
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6


# SGE_ID and JOD_ID
sge_task_id=0
if os.environ.has_key("SGE_TASK_ID"):
    # sge_task_id - 1 is to allow for easy indexing of lists...
    sge_task_id = long(os.environ["SGE_TASK_ID"]) - 1

job_id=0
if os.environ.has_key("JOB_ID"):
    job_id=long(os.environ["JOB_ID"])

print('Job id:', job_id)
print('Task id:', sge_task_id)

############################
# Establish SSH Connection #
############################
s = socket.socket()
try:
    s.connect(('localhost', 9002))
    print('Port forwarding to guybrush already exists')

except:
    print('Port forwarding to guybrush needs to be established')
    ssh_open_arg = ['ssh',
               '-4',
               '-M',
               '-S',
               'gurobi-socket',
               '-fnNT',
               '-L',
               '9002:guybrush-pi.compbio.ucsf.edu:41954',
               'sous.compbio.ucsf.edu']
    setup_port_forwarding = subprocess.Popen(ssh_open_arg)
    setup_port_forwarding.wait()

    s.connect(('localhost', 9002))
    print('Connection Established')

s.close()

#######################
# START OF GUROBI JOB #
#######################

time_start = roundTime()
print('Starting time:', time_start)

#################
# Set variables #
#################

print 'Starting Gurobi...'
from gurobipy import *

struct_id = int(sge_task_id + 1)
print('Working with conformer {0}'.format(struct_id))

#-- Check for command line arguements --#

# Get motif count to solve for from command line
motif_count = int(sys.argv[1])
print('Solving for {0} residue binding motifs (including ligand)'.format(motif_count))

extra_option = sys.argv[2] if len(sys.argv) > 2 else False

include = extra_option == 'include'
exclude = extra_option == 'exclude'

# If second argument is passed, it needs to be a text file of previous matches
# e.g. REN_0025-1_207_294_764_783

# Use residues from previous solutions in next iteration of solutions
if include:
    previous_match_list = [a for a in open(sys.argv[3], 'r')]
    print(previous_match_list)

    current_match_iteration = previous_match_list[sge_task_id].strip()
    current_match_total_split = re.split('-|_|', current_match_iteration)
    current_match_dash_split = re.split('-', current_match_iteration)

    print('Starting match iteration with {0}'.format(current_match_iteration))

    # struct_id is now whatever
    struct_id = int(current_match_total_split[1])
    print('struct_id set to {0}'.format(struct_id))

    # Get new residue constraints
    residue_constraints_string = current_match_dash_split[1]
    residue_constraints = [int(a) for a in residue_constraints_string.split('_')][1:]
    print('Adding residue constraints to optimization: {0}'.format(residue_constraints))

# Only use a subset of residues from fuzzball to solve for solutions
elif exclude:
    print('Only using a subset of fuzzball residues to generate solutions')
    include_position_resnumes = [int(a) for a in open(sys.argv[3], 'r')]

    print('Using resnums: {0}'.format(include_position_resnumes))

#-- Select pairwise scores for structid aka conformer --#

# Generate gurobi input table/csv
connection = sqlite3.connect('../two_body_terms.db')
cursor = connection.cursor()

residue_table = pd.read_sql_query("SELECT * from residues where struct_id = {0}".format(struct_id), connection)

# Upweight hbond_sc score term with ligand to bias towards solutions where hbond donor/acceptors on ligand are satisfied
# score_table = pd.read_sql_query(
#     """
#     SELECT struct_id, resNum1, resNum2,
#     CASE
#     WHEN sum(score_value) >= 0 THEN 1
#     WHEN sum(score_value) >= -0.1 THEN 0
#     ELSE sum(score_value)
#     END
#     as score_total from
#     (SELECT struct_id, resNum1, resNum2, score_type_name,
#     CASE
#     WHEN (score_type_name = 'hbond_sc' and resNum1 = 1) THEN 3 * score_value
#     ELSE score_value
#     END as score_value
#     from relevant_2b_scores where struct_id={0})
#     group by struct_id, resNum1, resNum2;
#     """.format(struct_id), connection)

score_table = pd.read_sql_query(
    """
    SELECT struct_id, resNum1, resNum2,
    CASE
    WHEN sum(score_value) >= 0.1 THEN 1
    ELSE sum(score_value)
    END
    as score_total from relevant_2b_scores where struct_id={0} group by struct_id, resNum1, resNum2;
    """.format(struct_id), connection)

#######################
# Solve for solutions #
#######################

# Set up model
residue_interactions = Model("residue_interactions")

# Get all possible residue pairs
MIP_var_dict = {}

# Add ligand
MIP_var_dict[float(1.0)] = residue_interactions.addVar(vtype=GRB.BINARY, name=str(1))

# Only add residues to Model if ligand-residue interaction energy is less than X
ligand_residue_scores = score_table.groupby(['struct_id', 'resNum1']).get_group((struct_id, 1))
for index, row in ligand_residue_scores.iterrows():
    if row['score_total'] < -0.5:
        MIP_var_dict[row['resNum2']] = residue_interactions.addVar(vtype=GRB.BINARY, name=str(row['resNum2']))

# List of residue indcies used in Model
MIP_residx_list = list(MIP_var_dict.keys())
print(ligand_residue_scores)
print(MIP_residx_list)
print(len(MIP_residx_list))

# Set up dict with pairwise scores
score_dict = {}
for index, row in score_table.iterrows():
    if all([row['resNum1'] in MIP_residx_list, row['resNum2'] in MIP_residx_list]):
        score_dict[(row['resNum1'], row['resNum2'])] = row['score_total']

#-- Set objective function --#
two_body_interactions = [MIP_var_dict[key[0]] * MIP_var_dict[key[1]] * value for key, value in score_dict.items()]
residue_interactions.setObjective(quicksum(two_body_interactions), GRB.MINIMIZE)

#-- Add constraints --#

# Always include ligand (residue 1)
residue_interactions.addConstr(MIP_var_dict[float(1)] == 1)

# Include residue constraints if sys.argv[2] exists
if include:
    for residue_num in residue_constraints:
        residue_interactions.addConstr(MIP_var_dict[float(residue_num)] == 1)

# Exclude residues not specified in posfile if posfile is provided
elif exclude:
    for resnum in MIP_residx_list:
        if resnum not in include_position_resnumes:
            residue_interactions.addConstr(MIP_var_dict[float(resnum)] == 0)

# Number of residues in a binding motif (includes ligand)
# NOTE: == for nucleation, <= for iterations!!!
residue_interactions.addConstr(quicksum(var for var in MIP_var_dict.values()) == motif_count)

# todo: update this to use score_dict
# Residues cannot be a solution if two-body interaction energy is above X
current_struct_scores = score_table.groupby(['struct_id']).get_group(struct_id)
for index, row in current_struct_scores.iterrows():
    if row['score_total'] > 0 and all([row['resNum1'] in MIP_residx_list, row['resNum2'] in MIP_residx_list]):
        residue_interactions.addConstr(MIP_var_dict[int(row['resNum1'])] + MIP_var_dict[int(row['resNum2'])] <= 1)

# --- SPECIAL CASES --- #

# Require one TRP in solution
# trp_residue_rows = residue_table[residue_table['res_type'] == 'TRP']
# print(trp_residue_rows)
# trp_resnum_list = [row['resNum'] for row in trp_residue_rows]
# residue_interactions.addConstr(quicksum(MIP_var_dict[res] for res in trp_resnum_list) == 1)

# Set Parameters
residue_interactions.Params.PoolSolutions = 50000
residue_interactions.Params.PoolGap = 0.5
residue_interactions.Params.PoolSearchMode = 2
residue_interactions.Params.Threads = 24

# May or may not be necessary...
# residue_interactions.Params.MIPFocus = 3
# residue_interactions.Params.Presolve = 2
# residue_interactions.Params.Heuristics = 0.1

# --- OPTIMIZE --- #
residue_interactions.optimize()

#-- Retrieve solutions --#

# Get residue-index mapping for current conformer
index_mapping = pd.read_sql_query("SELECT * from residue_index_mapping WHERE struct_id = {0}".format(struct_id), connection, index_col='residue_index')

# Set variables
results_list = []
compound_id = residue_table.loc[residue_table["resNum"] == 1]['name3'][0]
residues_in_motif = motif_count - 1

for i in range(residue_interactions.SolCount):
    residue_interactions.setParam(GRB.Param.SolutionNumber, i)

    res_index_tuple = [int(float(value.VarName)) for index, value in enumerate(residue_interactions.getVars()) if int(value.Xn) == 1]
    print('Residue Indicies: {}'.format(res_index_tuple))

    # Janky method to get values for non-ideal solutions since I can't get Model.PoolObjVal to work...
    solution_residue_pairs = [a for a in itertools.combinations(res_index_tuple, 2)]
    non_ideal_solution = sum(value for key, value in score_dict.items() if key in solution_residue_pairs)
    print('Non-Ideal Obj: {}'.format(non_ideal_solution))

    # todo: reorder lists so previous match residue indicies are at the end

    if include:
        new_residues = sorted(list(set(res_index_tuple) - set(residue_constraints)))
        res_index_tuple = new_residues + residue_constraints

    results_list.append({'Residue_indicies': res_index_tuple,
                         'Obj_score': non_ideal_solution,
                         'Conformer': '{}_{:0>4}'.format(compound_id, struct_id)})

df = pd.DataFrame(results_list)
if extra_option == 'include':
    df.to_csv('Gurobi_results-{0}-{1}_residue-conformer_{2}.csv'.format(compound_id, residues_in_motif, current_match_iteration))
else:
    df.to_csv('Gurobi_results-{0}-{1}_residue-conformer_{2}.csv'.format(compound_id, residues_in_motif, struct_id))

#####################
# END OF GUROBI JOB #
#####################

time_end = roundTime()
print('Ending time:', time_end)
print("Elapsed time:", time_end-time_start)

# Calculate RAM usage
qstat_p = subprocess.Popen(['/usr/local/sge/bin/linux-x64/qstat', '-j', '%d' % job_id], stdout=subprocess.PIPE)
out, err = qstat_p.communicate()

for line in out.split(os.linesep):
    m = re.match('(?:usage\s+%d[:]\s+.*?)(?:maxvmem[=])(\d+[.]\d+)([a-zA-Z]+)(?:.*?)' % sge_task_id, line)
    if m:
        ram_usage = float(m.group(1))
        ram_usage_type = m.group(2)
        print('Max virtual memory usage: %.1f%s' % (ram_usage, ram_usage_type))

#########################
# Close port forwarding #
#########################

ssh_close_arg = ['ssh',
                 '-S',
                 'gurobi_socket',
                 '-O',
                 'exit',
                 'sous.compbio.ucsf.edu']

close_port_forwarding = subprocess.Popen(ssh_close_arg)
close_port_forwarding.wait()