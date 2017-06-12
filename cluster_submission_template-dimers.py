#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-1624
#$ -l arch=linux-x64
#$ -l netapp=2G,scratch=1G

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

time_start = roundTime()

# SGE_ID and JOD_ID
sge_task_id=0
if os.environ.has_key("SGE_TASK_ID"):
    # sge_task_id - 1 is to allow for easy indexing of lists...
    sge_task_id = long(os.environ["SGE_TASK_ID"]) - 1

job_id=0
if os.environ.has_key("JOB_ID"):
    job_id=long(os.environ["JOB_ID"])

print('Starting time:', time_start)
print('Job id:', job_id)
print('Task id:', sge_task_id)

###########################
# Variable and Path Setup #
###########################

# Get Target Compound Directory as argv
target_compound_code = sys.argv[1].upper()

# Define relevant paths
bsff_path = os.path.join('/netapp', 'home', 'james.lucas', 'BindingSitesFromFragments')
scaffold_path = os.path.join(bsff_path, 'Dimer_Scaffolds')
target_compound_path = os.path.join(bsff_path, 'Compounds', target_compound_code)
constraint_file_path = os.path.join(target_compound_path, 'Pulled_Constraints')

# Get number of constraint files for target compound
number_of_cst_files = len(os.listdir(constraint_file_path))

# Get scaffold file name
# Modulus of SGE_ID by # of constraint files
scaffold_file = open(os.path.join(scaffold_path, 'filtered_heterodimers_biological_units_100.txt'), 'r')
scaffold_file_list = [file_name for file_name in scaffold_file]

current_scaffold = scaffold_file_list[sge_task_id]

current_scaffold_path = os.path.join(scaffold_path, 'cleaned_heterodimers_all_biological_units', current_scaffold)

# Get posfile and gridlig files for current scaffold
posfile_name = current_scaffold[:-3] + '.pos'
posfile_path = os.path.join(scaffold_path, 'posfiles', posfile_name)

gridlig_name = current_scaffold[:-3] + '.gridlig'
gridlig_path = os.path.join(scaffold_path, 'gridligs', gridlig_name)

# Get params file path
params_name = target_compound_code + '.params'
params_path = os.path.join(target_compound_path, params_name)

# Output path
output_path = os.path.join(target_compound_path, 'matches')

###########################
#   Construct Arguement   #
###########################

for cst_file in os.listdir(constraint_file_path):
    arg = ['/netapp/home/james.lucas/Rosetta/main/source/bin/match.linuxgccrelease',
           '-database',
           '/netapp/home/james.lucas/Rosetta/main/database',
           '-s',
           current_scaffold_path,
           '-match::lig_name',
           '{LIGAND THREE-LETTER CODE}',
           '-match::grid_boundary',
           gridlig_path,
           '-match::scaffold_active_site_residues',
           posfile_path,
           '-match::geometric_constraint_file',
           os.path.join(constraint_file_path, cst_file),
           '-extra_res_fa',
           params_path,
           '-output_matches_per_group',
           '1',
           '-match:consolidate_matches',
           '-ex1',
           '-ex2',
           '-extrachi_cutoff',
           '0',
           '-use_input_sc',
           '-euclid_bin_size',
           '1.5',
           '-euler_bin_size',
           '15',
           '-bump_tolerance',
           '0.5',
           '-out::path',
           output_path,
           '-match:output_format',
           'PDB']

    outfile_path = os.path.join(target_compound_path, 'stdout', 'rosetta.out')
    rosetta_outfile = open(outfile_path, 'w')

    rosetta_process = subprocess.Popen(arg, stdout=rosetta_outfile, cwd=os.getcwd())
    return_code = rosetta_process.wait()
    print('Task return code:', return_code, '\n')
    rosetta_outfile.close()

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

error_out = sys.argv[0] + '.e' + str(job_id) + '.' + str(sge_task_id)
output_out = sys.argv[0] + '.o' + str(job_id) + '.' + str(sge_task_id)

try:
    shutil.move(error_out , os.path.join(target_compound_path, 'stdout'))
    shutil.move(output_out ,  os.path.join(target_compound_path, 'stdout'))
except:
    print('No error or out file!')