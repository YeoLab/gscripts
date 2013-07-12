#!/bin/sh
#PBS -q glean
#PBS -N EVENT_TYPE
#PBS -l nodes=16:ppn:2
#PBS -l walltime=72:00:00
#PBS -o EVENT_TYPE_index/pbs.out
#PBS -e EVENT_TYPE_index/pbs.err
#PBS -V
#PBS -M olga.botvinnik@gmail.com
#PBS -m abe

# This is a template file for submitting gff indexing jobs
# to the TSCC PBS cluster. It is intended to be used by first replacing the text
# EVENT_TYPE with abbreviation the relevant event:
# - Skipped exons (SE)
# - Alternative 3’/5’ splice sites (A3SS, A5SS)
# - Mutually exclusive exons (MXE)
# - Tandem 3’ UTRs (TandemUTR)
# - Retained introns (RI)
# - Alternative first exons (AFE)
# - Alternative last exons (ALE)

# if not all the chromosomes are there in the folder, remove the folder
# because otherwise MISO will detect a non empty folder and abort.
# Check for existence of Chr9 because it's one of the last ones that gets indexed.
if [ -d

DIR=/home/obotvinnik/genome/miso_annotations/hg19
OUT_DIR=$DIR/EVENT_TYPE_index/

python /home/yeo-lab/software/bin/index_gff.py --index $DIR/EVENT_TYPE.hg19.gff3 $DIR/EVENT_TYPE_index/