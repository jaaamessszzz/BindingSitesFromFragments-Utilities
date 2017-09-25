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

# DEBUGGING
get_env = subprocess.Popen(['env'])
get_env.wait()

from gurobipy import *

struct_id = int(sge_task_id + 1)

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

