#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -pe smp 30
#$ -R yes
#$ -l h_rt=240:00:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=5G
#$ -l netapp=5G,scratch=5G

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

##########################
# Start submitting tasks #
##########################

time_start = roundTime()
print('Starting time:', time_start)

arg = ['scl',
       'enable',
       'python27',
       '\"python /netapp/home/james.lucas/BindingSitesFromFragments/BindingSitesFromFragments-Utilities/Gurobi/gurobi_cluster-job.py {0}\"'.format(sge_task_id + 1),
       ]

print(' '.join(arg))

outfile_path = os.path.join('stdout', 'gurobi_out-{0}.out'.format(sge_task_id))
gurobi_outfile = open(outfile_path, 'w')
gurobi_process = subprocess.Popen(arg, stdout=gurobi_outfile, cwd=os.getcwd())
return_code = gurobi_process.wait()

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

