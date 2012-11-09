#!/bin/sh
#$ -N score 
#$ -S /bin/sh
#$ -o score.out
#$ -e score.log
#$ -l bigmem
#$ -l h_vmem=40g
#$ -j n
#$ -cwd
#$ -V
#$ -p 0

for i in {1...100}
do
  echo $i
  python /nas3/gpratt/projects/structure/RNAse/get_regions.py /nas3/gpratt/projects/structure/RNAse/hela_xlink_noprotein.plus.strscore.bw /nas3/gpratt/projects/structure/RNAse/hela_xlink_noprotein.minus.strscore.bw shuffled_regions_all${i}_v3.bed /nas3/gpratt/ens/hg19.genome > shuffled_regions_all${i}_scored_v3.bed 

done