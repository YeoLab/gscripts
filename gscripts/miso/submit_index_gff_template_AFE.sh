#!/bin/sh
#PBS -q glean
#PBS -N AFE
#PBS -l nodes=16:ppn:2
#PBS -l walltime=72:00:00
#PBS -o /home/obotvinnik/genome/miso_annotations/hg19/AFE_index/pbs.out
#PBS -e /home/obotvinnik/genome/miso_annotations/hg19/AFE_index/pbs.err
#PBS -V
#PBS -M olga.botvinnik@gmail.com
#PBS -m abe

# This is a template file for submitting gff indexing jobs
# to the TSCC PBS cluster. It is intended to be used by first replacing the text
# AFE with abbreviation the relevant event:
# - Skipped exons (SE)
# - Alternative 3’/5’ splice sites (A3SS, A5SS)
# - Mutually exclusive exons (MXE)
# - Tandem 3’ UTRs (TandemUTR)
# - Retained introns (RI)
# - Alternative first exons (AFE)
# - Alternative last exons (ALE)
# and then submitted. SO for example,
# PREFIX=/home/obotvinnik/gscripts/gscripts/miso/submit_index_gff_template 
# AFE=AFE 
# SUBMIT_SH=$PREFIX\_$AFE.sh 
# sed "s:AFE:$AFE:g" $PREFIX.sh > $SUBMIT_SH 
# qsub $SUBMIT_SH
#
# Or, all one one line for copying:
# PREFIX=/home/obotvinnik/gscripts/gscripts/miso/submit_index_gff_template ; AFE=AFE ; SUBMIT_SH=$PREFIX\_$AFE.sh ; sed "s:AFE:$AFE:g" $PREFIX.sh > $SUBMIT_SH ; qsub $SUBMIT_SH

DIR=/home/obotvinnik/genome/miso_annotations/hg19
OUT_DIR=$DIR/AFE_index/

# if not all the chromosomes are there in the folder, remove the folder
# because otherwise MISO will detect a non empty folder and abort.
# Check for existence of Chr9 because it's one of the last ones that gets indexed.
if [ -d $OUT_DIR ] ; then
    if [ -d $OUT_DIR/chr9 ] ; then
	rm -rf $OUT_DIR
    fi
fi

python /home/yeo-lab/software/bin/index_gff.py --index $DIR/AFE.hg19.gff3 $DIR/AFE_index/