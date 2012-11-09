#!/bin/sh

#This works for mRNA mapped peaks only
DIR=/nas3/gpratt
for f in $@
do
  echo $f;
  annotatePeaks.pl $f hg19 -size 4000 -hist 10 -d $DIR/RBP_peaks/hnrnpa1 $DIR/RBP_peaks/hnrnpa2b1 $DIR/RBP_peaks/hnrnpf $DIR/RBP_peaks/hnrnph1 $DIR/RBP_peaks/hnrnpm $DIR/RBP_peaks/hnrnpu > $f.density;
done