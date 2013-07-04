#!/bin/sh

for f in /nas3/gpratt/RBP_peaks/*liftover_hg19.bed
do
  NAME=${f%*.*}
  ./ensMapper.sh $f > $NAME.ens.bed
  intersectBed -a $NAME.ens.bed -b /nas3/gpratt/ens/ensembl_genes_ens_cds_start_end_padded100.bed > $NAME.ens.padded100.bed
done