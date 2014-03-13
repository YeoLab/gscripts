#!/bin/sh
#PBS -N test_qtools_submitter
#PBS -o test_qtools_submitter.out
#PBS -e test_qtools_submitter.err
#PBS -V
#PBS -l walltime=0:01:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home-yeo

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
date
echo testing PBS

