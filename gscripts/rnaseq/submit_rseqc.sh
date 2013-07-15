#!/bin/sh
#PBS -q glean
#PBS -N rnaseqc
#PBS -l nodes=16:ppn=2
#PBS -l walltime=36:00:00
#PBS -o BASE_DIR/rseqc.out
#PBS -e /home/obotvinnik/genomes/miso_annotations/hg19/submit_index_gff_template_EVENT_TYPE.err
#PBS -V
#PBS -M olga.botvinnik@gmail.com
#PBS -m abe

DIR=BASE_DIR
RNASEQC=$DIR/rnaseqc/

java -jar /home/yeo-lab/software/rnaseqc/RNA-SeQC_v1.1.7.jar -o $RNASEQC_DIR