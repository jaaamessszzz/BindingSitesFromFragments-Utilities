#! /usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=00:10:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=1G
#$ -l netapp=1G,scratch=1G

import subprocess

# ssh port forwarding
ssh_open = ['ssh',
            '-4',
            '-M',
            '-S',
            'gurobi_socket',
            '-fnNT',
            '-L',
            '9002:guybrush-pi.compbio.ucsf.edu:41954',
            'sous.compbio.ucsf.edu'
            ]

ssh_open_p = subprocess.Popen(ssh_open)
ssh_open_p.wait()

# Gurobi shell
gurobi_shell = subprocess.Popen(['scl',
                                 'enable',
                                 'python27',
                                 'gurobi.sh'
                                 ]
                                )
gurobi_shell.wait()

# Gurobi python import
gurobi_python = subprocess.Popen(['scl',
                                  'enable',
                                  'python27',
                                  'python /netapp/home/james.lucas/BindingSitesFromFragments/BindingSitesFromFragments-Utilities/import_gurobi.py'
                                  ]
                                 )
gurobi_python.wait()

# close ssh
ssh_close = ['ssh',
             '-S',
             'gurobi_socket',
             '-O',
             'exit',
             'sous.compbio.ucsf.edu'
             ]
