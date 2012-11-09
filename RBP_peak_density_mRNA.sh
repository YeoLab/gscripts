#!/bin/sh

#This works for mRNA mapped peaks only
DIR=/nas3/gpratt
for f in $@
do
  echo $f;
  annotatePeaks.pl $f $DIR/ens/ensembl_genes.fasta -size 4000 -hist 10 -d $DIR/RBP_peaks/hnrnpa1.clip.293t.hg19.bowtie.ens $DIR/RBP_peaks/hnrnpa2b1.clip.293t.hg19.bowtie.ens $DIR/RBP_peaks/hnrnpf.clip.293t.hg19.bowtie.ens $DIR/RBP_peaks/hnrnph1.clip.293t.hg19.bowtie.ens $DIR/RBP_peaks/hnrnpm.clip.293t.hg19.bowtie.ens $DIR/RBP_peaks/hnrnpu.clip.293t.hg19.bowtie.ens > $f.density;
done