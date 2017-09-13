#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r yes
#$ -l h_rt=00:10:00
#$ -t 1
#$ -l arch=linux-x64
#$ -l mem_free=1G
#$ -l netapp=1G,scratch=1G

ssh -4 -M -S gurobi_socket -fnNT -L 9002:guybrush-pi.compbio.ucsf.edu:41954 sous

scl enable python27 'python /netapp/home/james.lucas/BindingSitesFromFragments/BindingSitesFromFragments-Utilities/import_gurobi.py'

ssh -S gurobi_socket -O exit sous