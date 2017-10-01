#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24
#$ -R yes
#$ -l h_rt=240:00:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=10G
#$ -l netapp=5G,scratch=1G
#$ -l q=lab.q

# Yes, I really need 240GB of RAM
scl enable python27 'python /netapp/home/james.lucas/BindingSitesFromFragments/BindingSitesFromFragments-Utilities/Gurobi/gurobi_cluster-complete_job.py'