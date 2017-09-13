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

try:
    gurobi_shell = subprocess.Popen(['scl',
                                     'enable',
                                     'python27',
                                     'gurobi.sh'
                                     ]
                                    )
    gurobi_shell.wait()

except Exception as e:
    print "Shell failed with exception: {0}".format(e)

try:
    gurobi_python = subprocess.Popen(['scl',
                                      'enable',
                                      'python27',
                                      'python /netapp/home/james.lucas/BindingSitesFromFragments/BindingSitesFromFragments-Utilities/import_gurobi.py'
                                      ]
                                     )
    gurobi_python.wait()

except Exception as e:
    print "Python Import faied with exception {0}".format(e)

ssh_close = ['ssh',
             '-S',
             'gurobi_socket',
             '-O',
             'exit',
             'sous.compbio.ucsf.edu'
             ]
