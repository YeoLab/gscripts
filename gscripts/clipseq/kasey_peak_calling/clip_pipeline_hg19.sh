#!/bin/sh
#
# ------------------
#$ -N FOO
#$ -S /bin/sh
##$ -o .out
##$ -e qsub.log
##$ -M khutt@ucsd.edu
#$ -m es
# ---------------------------
#
#$ -V
# Execute the job from  the  current  working  directory
#$ -cwd
#
# All resources are defined here
#
# Choose your queue
#
##$ -q <queuename>
#
# Job priority
#
#$ -p 0.1
#$ -l h_vmem=60G
# ---------------------------
#
# Put compilations here
#
# ---------------------------
#
# Execution
#

cutadapt -O 5 -b TCGTATGCCGTCTTCTGCTTG -o foo_out foo_in
perl /nas3/shuelga/gscripts/run_bowtie12.7_20seed_qual.pl hg19 foo_out
perl /nas3/shuelga/gscripts/run_bowtie12.7_20seed_qual_repeats.pl hg19 foo_out
perl /nas3/shuelga/gscripts/readgenomealignments2_reomverepeats.pl hg19 foo_out
perl /nas3/oootif/Solexa/scripts/YL_combine_reads2.pl foo_out foo_out.unique
##perl /nas3/khutt/DLC_sequencing/FET_clip/MCF7/allocate_reads_to_genes.pl hg19 foo.output
perl /nas3/khutt/scripts_peak_calling/allocate_reads_to_genes_hg19.pl hg19 foo_out
perl /nas3/oootif/Solexa/scripts/BED_combine_sort.pl hg19 foo_out.ingenes.BED foo_out.notingenes.BED > foo_out.all.BED
perl /nas3/khutt/wig_script/BED_to_wiggle_fast_hg18.pl foo_out.all.BED 255,0,0
perl /nas3/khutt/scripts_peak_calling/analyze_clip_data_complete hg19 foo_out.ingenes.BED -force -shift 0 -trim 0 -window_min 50 -window_max 150 -gap 20 -mRNA 0 -pval 0.0001 0
##perl /nas3/khutt/scripts_peak_calling/analyze_clip_data_complete mm9 TAF15_195.ingenes.BED -force -shift 0 -trim 1 -mRNA 0
