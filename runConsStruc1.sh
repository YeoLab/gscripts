#!/bin/sh 
#
# ------------------
#$ -N all_cons
#$ -S /bin/sh
#$ -o dise.out
#$ -e dise.log
#$ -l bigmem
#$ -j n
# ---------------------------
#
# Execute the job from  the  current  working  directory
#$ -cwd
#$ -V
# All resources are defined here
#
# Choose your queue
#
##$ -q <queuename>
#
# Job priority
#
#$ -p 0
#
# ---------------------------
#
# Put compilations here
#
# ---------------------------
#
# Execution
#
echo "hostname is:"
hostname
python /nas3/lovci/gscripts/conserved_structure.py --aligned_links /nas3/lovci/projects/structure/hg19/ultra_conserved_links.txt --all_links /nas3/gpratt/gscripts/shuffled_regions_all.bed > conserved_structure.all.shuffled.txt