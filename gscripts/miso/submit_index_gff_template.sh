#!/bin/sh
#PBS -q glean
#PBS -N EVENT_TYPE
#PBS -l nodes=16:ppn=2
#PBS -l walltime=36:00:00
#PBS -o /home/obotvinnik/genomes/miso_annotations/hg19/submit_index_gff_template_EVENT_TYPE.out
#PBS -e /home/obotvinnik/genomes/miso_annotations/hg19/submit_index_gff_template_EVENT_TYPE.err
#PBS -V
#PBS -M olga.botvinnik@gmail.com
#PBS -m abe

# Date: 07-12-2013
# Author: Olga Botvinnik
#
# This is a template file for submitting miso gff indexing jobs
# to the TSCC PBS cluster. It is intended to be used by first replacing the text
# EVENT_TYPE with abbreviation the relevant event:
# - Skipped exons (SE)
# - Alternative 3’/5’ splice sites (A3SS, A5SS)
# - Mutually exclusive exons (MXE)
# - Tandem 3’ UTRs (TandemUTR)
# - Retained introns (RI)
# - Alternative first exons (AFE)
# - Alternative last exons (ALE)
# and then submitted. SO for example,
# PREFIX=/home/obotvinnik/gscripts/gscripts/miso/submit_index_gff_template 
# EVENT_TYPE=AFE 
# SUBMIT_SH=$PREFIX\_$EVENT_TYPE.sh 
# sed "s:EVENT_TYPE:$EVENT_TYPE:g" $PREFIX.sh > $SUBMIT_SH 
# qsub $SUBMIT_SH
#
# Or, all one one line for copying:
# TEMPLATE_DIR=/home/obotvinnik/gscripts/gscripts/miso ; DIR=/home/obotvinnik/genomes/miso_annotations/hg19 ; PREFIX=submit_index_gff_template ; EVENT_TYPE=AFE ; TEMPLATE_SH=$TEMPLATE_DIR/$PREFIX.sh ; SUBMIT_SH=$DIR/$PREFIX\_$EVENT_TYPE.sh ; sed "s:EVENT_TYPE:$EVENT_TYPE:g" $TEMPLATE_SH > $SUBMIT_SH ; qsub $SUBMIT_SH
# to run on all event types:
# EVENT_TYPES="A3SS A5SS AFE ALE MXE SE TandemUTR RI" ; TEMPLATE_DIR=/home/obotvinnik/gscripts/gscripts/miso ; DIR=/home/obotvinnik/genomes/miso_annotations/hg19 ; PREFIX=submit_index_gff_template ; TEMPLATE_SH=$TEMPLATE_DIR/$PREFIX.sh ; for EVENT_TYPE in $EVENT_TYPES ; do SUBMIT_SH=$DIR/$PREFIX\_$EVENT_TYPE.sh ; sed "s:EVENT_TYPE:$EVENT_TYPE:g" $TEMPLATE_SH > $SUBMIT_SH ; qsub $SUBMIT_SH ; done

source /home/obotvinnik/.bashrc

DIR=/home/obotvinnik/genomes/miso_annotations/hg19
OUT_DIR=$DIR/EVENT_TYPE_index/

# if not all the chromosomes are there in the folder, remove the folder
# because otherwise MISO will detect a non empty folder and abort.
# Check for existence of Chr9 because it's one of the last ones that gets indexed.
if [ -d $OUT_DIR ] ; then
    if [ ! -d $OUT_DIR/chr9 ] ; then
	rm -rf $OUT_DIR
    fi
fi

python /home/yeo-lab/software/bin/index_gff.py --index $DIR/EVENT_TYPE.hg19.gff3 $DIR/EVENT_TYPE_index/