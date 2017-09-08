#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=00:10:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=1G
#$ -l netapp=1G,scratch=1G

import sys
import subprocess
import os

ssh_open = ['ssh',
            '-M',
            '-S',
            'my-ctrl-socket',
            '-fnNT',
            '-L',
            '7001:guybrush.ucsf.edu:41954',
            'sous.compbio.ucsf.edu'
            ]

ssh_open_p = subprocess.Popen(ssh_open)
ssh_open_p.wait()
gurobi_test = subprocess.Popen(['gurobi.sh'])
gurobi_test.wait()

ssh_close = ['ssh',
             '-S',
             'my-ctrl-socket',
             '-O',
             'exit',
             'sous.compbio.ucsf.edu'
             ]