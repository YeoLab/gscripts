#!/bin/sh
#$ -N shuffler 
#$ -S /bin/sh
#$ -o mfe.out
#$ -e mfe.log
#$ -l bigmem
#$ -l h_vmem=40g
#$ -j n
#$ -cwd
#$ -V
#$ -p 0

for i in {1..100}
do
  echo $i
  python /nas3/lovci/gscripts/conserved_structure.py --aligned_links /nas3/lovci/projects/structure/hg19/ultra_conserved_links.txt --all_links shuffled_regions_all${i}_v3.bed > conserved_regions_shuffled${i}_v3.bed
done