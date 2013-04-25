#!/bin/sh

#Takes any number of bam files and outputs them as bed files mapped to the ensomble genes
for f in $@
do
  FILE=${f%*.*}
  OUTPUTFILE=`echo $FILE | sed 's/hg19/ens/g'`
 
  intersectBed -bed -wo -abam $FILE.bam -b /nas3/gpratt/ens/ensembl_genes_exons.bed > $FILE.bed
 
  cat $FILE.bed | awk 'BEGIN {OFS="\t"} { if($10=="+") {print $11, ($12 + ($2-$8)), ($12 + ($3 - $8))} else {print $11, ($12 + ($9-$3)), ($12 + ($9-$2))}}' | awk 'BEGIN {OFS="\t"} {if($2 < 0) {$2 = 0}; if($3 < 0) {$3 = 1}; print $0}'> $OUTPUTFILE.bed
done 