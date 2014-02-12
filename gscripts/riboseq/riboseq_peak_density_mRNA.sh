#!/bin/sh

for f in $@
do
  echo $f;
  annotatePeaks.pl $f /nas3/gpratt/ens/ensembl_genes.fasta -size 4000 -hist 10 -d /nas3/gpratt/HEK293/GSM782786_NAMF-mapped.ens /nas3/gpratt/HEK293/GSM782788_NACF-mapped.ens > $f.density;
done