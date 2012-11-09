#!/bin/sh

#Takes any number of bed files and outputs them as bed files mapped to the ensomble genes

#need to make sure files have exactly 6 columns
f=$1
FILE=${f%*.*}
EXONS=$2

#/nas3/gpratt/ens/ensembl_genes_exons.bed need to fix this in human later
cat $f | awk 'BEGIN {OFS="\t"} {print $1,$2,$2,NR,$1,$1,$1,$1}' | cut -f 1,2,3,4,5,6 > $f.start.cut

cat $f | awk 'BEGIN {OFS="\t"} {print $1,$3,$3,NR,$1,$1,$1,$1}' | cut -f 1,2,3,4,5,6 > $f.stop.cut

intersectBed -wo -a $f.start.cut -b $2 > ${f}_intersected.start.bed
intersectBed -wo -a $f.stop.cut  -b $2 > ${f}_intersected.stop.bed


#do a sanity check to make sure everything matches up
startLen=`wc -l ${f}_intersected.start.bed | cut -d " " -f 1`
endLen=`wc -l ${f}_intersected.stop.bed | cut -d " " -f 1`

cat ${f}_intersected.start.bed | awk 'BEGIN {OFS="\t"} { if($10=="+") {print $11, ($12 + ($2-$8)), ($12 + ($3 - $8)), $4, $10} else {print $11, ($12 + ($9-$3)), ($12 + ($9-$2)), $4, $10}}' | awk 'BEGIN {OFS="\t"} {if($2 < 0) {$2 = 0}; if($3 < 0) {$3 = 1}; print $0}' | sort -n -k 4 > start.tmp.bed

cat ${f}_intersected.stop.bed | awk 'BEGIN {OFS="\t"} { if($10=="+") {print $11, ($12 + ($2-$8)), ($12 + ($3 - $8)), $4, $10} else {print $11, ($12 + ($9-$3)), ($12 + ($9-$2)), $4, $10}}' | awk 'BEGIN {OFS="\t"} {if($2 < 0) {$2 = 0}; if($3 < 0) {$3 = 1}; print $0}'  | sort -n -k 4 > stop.tmp.bed

#I'm loosing a lot of reads over here, need to fix this eventually...
python output.py start.tmp.bed stop.tmp.bed 
#rm ${f}_intersected.bed
#rm $f.cut
