#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-20000
#$ -l arch=linux-x64
#$ -l netapp=6G,scratch=1G

# *** qsub from BindingSitesFromFragments root directory *** #
# Make sure you set task number above to be correct!!!

# TEP cst files - 655
# 100% identity dimer scaffolds - 1624

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

target_compound_code = sys.argv[1]
block_size = int(sys.argv[2])

##########################
# Start submitting tasks #
##########################

matcher_arg_json = json.load(open('matcher_argument_list.json', 'r'))
current_arg_block = matcher_arg_json[(sge_task_id * block_size):((sge_task_id + 1) * block_size)]

bsff_path = os.path.join('/netapp', 'home', 'james.lucas', 'BindingSitesFromFragments')
target_compound_path = os.path.join(bsff_path, 'Compounds', target_compound_code)

for block in current_arg_block:
    time_start = roundTime()
    print('Starting time:', time_start)

    arg = ['/netapp/home/james.lucas/Rosetta/main/source/bin/match.linuxgccrelease',
           '-database',
           '/netapp/home/james.lucas/Rosetta/main/database',
           '-s',
           block[0],
           '-match::lig_name',
           block[1],
           '-match::grid_boundary',
           block[2],
           '-match::scaffold_active_site_residues',
           block[3],
           '-match::geometric_constraint_file',
           block[4],
           '-extra_res_fa',
           block[5],
           '-output_matches_per_group',
           '1',
           '-match:consolidate_matches',
           '-ex1',
           '-ex2',
           '-extrachi_cutoff',
           '0',
           '-use_input_sc',
           '-euclid_bin_size',
           '2', # Roland : 1.5
           '-euler_bin_size',
           '20', # Roland: 15
           '-bump_tolerance',
           '0.5',
           '-out::path',
           block[6],
           '-match:output_format',
           'PDB',
           '-out:file:scorefile', # match scores
           '-match_grouper',
           'SameSequenceGrouper', # Two matches belong in the same group if their hits come from the same amino acds at the same scaffold build positions
           '-mute',
           'protocols.idealize']

    print(' '.join(arg))

    scaffold = os.path.basename(os.path.normpath(block[0])).split('.')[0]
    constraint = os.path.basename(os.path.normpath(block[4])).split('.')[0]

    outfile_path = os.path.join(os.path.split(block[6])[0], 'stdout', '{0}-{1}.out'.format(scaffold, constraint))
    rosetta_outfile = open(outfile_path, 'w')
    rosetta_process = subprocess.Popen(arg, stdout=rosetta_outfile, cwd=os.getcwd())
    return_code = rosetta_process.wait()

    print('Task return code:', return_code, '\n')

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

error_out = '{0}.e{1}.{2}'.format(sys.argv[0], str(job_id), str(sge_task_id))
output_out = '{0}.o{1}.{2}'.format(sys.argv[0], str(job_id), str(sge_task_id))

print(error_out)
print(output_out)

try:
    shutil.move(error_out, os.path.join(target_compound_path, 'stdout'))
    shutil.move(output_out, os.path.join(target_compound_path, 'stdout'))
except:
    print('No error or out file!')