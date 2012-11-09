#!/bin/sh

for f in $@
do
  echo $f;
  annotatePeaks.pl $f hg19 -size 4000 -hist 10 -d /nas3/gpratt/HEK293/GSM782786_NAMF-mapped.hg19 /nas3/gpratt/HEK293/GSM782788_NACF-mapped.hg19 > $f.density;
done