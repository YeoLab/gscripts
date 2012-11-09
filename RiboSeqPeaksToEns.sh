#!/bin/sh

for f in /nas3/gpratt/HEK293/*.hg19.bam.filtered05
do
  OUTFILE=`echo $f | sed 's/hg19/ens/g'`
  echo $OUTFILE
  ./ensMapper.sh $f > $OUTFILE
  intersectBed -a ${OUTFILE} -b /nas3/gpratt/ens/ensembl_genes_ens_cds_start_end_padded100.bed > $OUTFILE.padded100.bed
done