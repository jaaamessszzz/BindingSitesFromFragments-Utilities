#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 28
#$ -R yes
#$ -l h_rt=240:00:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=5G
#$ -l netapp=5G,scratch=5G

scl enable python27 'python /netapp/home/james.lucas/BindingSitesFromFragments/BindingSitesFromFragments-Utilities/Gurobi/gurobi_cluster-complete_job.py'