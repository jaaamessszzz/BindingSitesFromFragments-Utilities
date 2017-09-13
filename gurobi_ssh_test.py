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
            'gurobi_socket',
            '-fnNT',
            '-L',
            '9002:guybrush.ucsf.edu:41954',
            'sous.compbio.ucsf.edu'
            ]

ssh_open_p = subprocess.Popen(ssh_open)
ssh_open_p.wait()
gurobi_test = subprocess.Popen(['scl',
                                'enable',
                                'python27',
                                '\'gurobi.sh\''
                                ]
                               )
gurobi_test.wait()

ssh_close = ['ssh',
             '-S',
             'gurobi_socket',
             '-O',
             'exit',
             'sous.compbio.ucsf.edu'
             ]
