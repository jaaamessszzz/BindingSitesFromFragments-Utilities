#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -R yes
#$ -j yes
#$ -l h_rt=240:00:00
#$ -t 1-10
#$ -l arch=linux-x64
#$ -l mem_free=5G
#$ -l netapp=1G,scratch=1G

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

def determine_matched_residue_positions(match_pdb_path):
    """
    Parse the filename of the match PDB to determine IDs and positions of match residues
    :return: 
    """
    positions_block = os.path.basename(os.path.normpath(match_pdb_path)).split('_')[2]
    resnames = [a for a in re.split("[0-9]*", positions_block) if a]
    resnums = [int(a) for a in re.split("[a-zA-Z]*", positions_block) if a]

    return [(a, b) for a, b in zip(resnames, resnums)]

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

input_pdb = sys.argv[1]
params_file_path = sys.argv[2]
design_json_path = sys.argv[3]
designs_per_task = sys.argv[4]

input_pdb_base = os.path.basename(os.path.normpath(input_pdb))
design_json = json.load(open(design_json_path, 'r'))

##########################
# Start submitting tasks #
##########################

time_start = roundTime()
print('Starting time:', time_start)

print(input_pdb_base)
print(determine_matched_residue_positions(input_pdb_base))
print(params_file_path)
print(design_json)
print(','.join([str(a) for a in design_json['design_residue_list']]))

for design in range(1, designs_per_task+1):
    arg = ['/netapp/home/james.lucas/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease',
           '-database',
           '/netapp/home/james.lucas/Rosetta/main/database',
           '-s',
           input_pdb,
           '-extra_res_fa',
           params_file_path,
           '-parser:protocol',
           'Design_Template.xml',
           # '-use_input_sc',
           # '-flip_HNQ',
           # '-no_optH',
           # 'false',
           '-holes:dalphaball', # Required for RosettaHoles Filter
           '/netapp/home/james.lucas/Rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc', # Full path to DAlphaBall.gcc on chef
           '-out:prefix',
           '{0}-'.format(sge_task_id + 1),
            '-out:suffix',
           '-{0}'.format(design),
           '-parser:script_vars',
           'motif_residues={0}'.format(','.join([str(resnum[1]) for resnum in determine_matched_residue_positions(input_pdb_base)])),
           'design_positions={0}'.format(','.join([str(a) for a in design_json['design_residue_list']])),
           'design_xml={0}'.format('BackrubEnsemble-Design.xml')
           ]

    print(' '.join(arg))

    rosetta_process = subprocess.Popen(arg, cwd=os.getcwd())
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

