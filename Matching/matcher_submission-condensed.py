#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-1000
#$ -j y
#$ -l arch=linux-x64
#$ -l mem_free=15G
#$ -l netapp=5G,scratch=1G

# Make sure you set task number above to be correct!!!

import socket
import sys
import subprocess
import os
from datetime import *
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

########################
# Positional arguments #
########################

block_size = int(sys.argv[1])

##########################
# Start submitting tasks #
##########################

matcher_arg_json = json.load(open('matcher_argument_list.json', 'r'))
current_arg_block = matcher_arg_json['args'][(sge_task_id * block_size):((sge_task_id + 1) * block_size)]

# Make directory for each job and switch into it... this is to prevent multiple simultaneous writes to matcher_score.sc
task_id_dir = 'taskid-{0}'.format(sge_task_id + 1)
os.mkdir(task_id_dir)
os.chdir(task_id_dir)

constraint_json_name = '{0}-constraint_blocks.json'.format(matcher_arg_json["target_compound"])
print(constraint_json_name)

# --- Prepare directory for task-specific constraint files --- #

if os.path.exists(constraint_json_name):
    print('{0} found, preparing directory for constraint files...'.format(constraint_json_name))
    os.mkdir(os.path.join(matcher_arg_json['cst_dir'], task_id_dir))

for block in current_arg_block:
    time_start = roundTime()
    print('Starting time:', time_start)

    # For reference...
    # block = [conformer, constraint, current_scaffold, posfile_name, params_name] + [gridlig]

    # --- Generate constraint file from constraint json and dump into matcher_arg_json['cst_dir'] --- #

    cst_file_path = os.path.join(matcher_arg_json['cst_dir'], block[1])

    if os.path.exists(constraint_json_name):
        print('{0} found, generating constraint file...'.format(constraint_json_name))
        with open(constraint_json_name, 'r') as cst_json:
            # Get required blocks from constraint file name
            cst_json_json = json.load(cst_json)
            relevant_constraint_blocks = re.split('-|\.', block[1])[1].split('_')[1:]

            cst_file_path = os.path.join(matcher_arg_json['cst_dir'], task_id_dir, block[1])

            with open(cst_file_path, 'w') as cst_file:
                for cst_block_index in relevant_constraint_blocks:
                    cst_file.write('CST::BEGIN\n')
                    cst_file.write(cst_json_json[block[0]][cst_block_index])
                    cst_file.write('\n')
                    cst_file.write('CST::END\n')

    arg = ['/netapp/home/james.lucas/Rosetta/main/source/bin/match.linuxgccrelease',
           '-database',
           '/netapp/home/james.lucas/Rosetta/main/database',
           '-s',
           os.path.join(matcher_arg_json['scaffold_pdb_dir'], block[2]),
           '-match::lig_name',
           matcher_arg_json['target_compound'],
           '-match::scaffold_active_site_residues',
           os.path.join(matcher_arg_json['posfile_dir'], block[3]),
           '-match::geometric_constraint_file',
           cst_file_path,
           '-extra_res_fa',
           os.path.join(matcher_arg_json['params_dir'], block[4]),
           '-output_matches_per_group',
           '1',
           '-match:consolidate_matches',
           '-ex1',
           '-ex2',
           '-extrachi_cutoff',
           '0',
           '-use_input_sc',
           '-euclid_bin_size',
           '1', # Roland : 1.5
           '-euler_bin_size',
           '10', # Roland: 15
           '-bump_tolerance',
           '0.5',
           '-out::path',
           matcher_arg_json['output_path'],
           '-match:output_format',
           'PDB',
           '-out:file:scorefile', # match scores
           '-match_grouper',
           'SameSequenceGrouper', # Two matches belong in the same group if their hits come from the same amino acids at the same scaffold build positions
           '-mute',
           'protocols.idealize']

    if len(block) == 7:
        arg.append('-match::grid_boundary')
        arg.append(os.path.join(matcher_arg_json['gridlig_path'], block[5]))

    print(' '.join(arg))

    rosetta_process = subprocess.Popen(arg, cwd=os.getcwd())
    return_code = rosetta_process.wait()

    print('Task return code:', return_code, '\n')

    # Delete constraint file after it is done being used
    os.remove(cst_file_path)

    time_end = roundTime()
    print('Ending time:', time_end)
    print("Elapsed time:", time_end-time_start)

    # Calculate RAM usage
    # qstat_p = subprocess.Popen(['/usr/local/sge/bin/linux-x64/qstat', '-j', '%d' % job_id], stdout=subprocess.PIPE)
    # out, err = qstat_p.communicate()

    # for line in out.split(os.linesep):
    #     m = re.match('(?:usage\s+%d[:]\s+.*?)(?:maxvmem[=])(\d+[.]\d+)([a-zA-Z]+)(?:.*?)' % sge_task_id, line)
    #     if m:
    #         ram_usage = float(m.group(1))
    #         ram_usage_type = m.group(2)
    #         print('Max virtual memory usage: %.1f%s' % (ram_usage, ram_usage_type))