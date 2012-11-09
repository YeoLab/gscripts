#!/bin/sh
#$ -N shuffler 
#$ -S /bin/sh
#$ -o shuff.out
#$ -e shuff.log
#$ -l bigmem
#$ -l h_vmem=40g
#$ -j n
#$ -cwd
#$ -V
#$ -p 0

for i in {1..100}
do
  echo $i
  python shuffle_bridges.py > shuffled_regions_all${i}_v3.bed
done