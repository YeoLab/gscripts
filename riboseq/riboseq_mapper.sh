#!/bin/bash

f=$1
NAME=${f%*.*}
fastq_collapser $NAME.fastq > $NAME.collapsed.fastq
fastq_quality_trimmer -Q33 -t 13 -l 20 -i $NAME.collapsed.fastq -o $NAME.trimmed.fastq
fastx_clipper -Q33 -a CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTGAA -i $NAME.trimmed.fastq -o $NAME.clipped.fastq
bowtie -S -p 4 -l 26 --un $NAME.not_aligned_rRNA.fasta /nas3/gpratt/bowtie/rRNA_transcripts $NAME.clipped.fastq $NAME.rRNA_aligned

bowtie -m 15 --best -p 4 --strata -S -l 26 --un $NAME.not_aligned_rRNA_all_ref.fasta all_ref $NAME.not_aligned_rRNA.fasta $NAME.rep_aligned.sam

bowtie -m 15 --best -p 4 --strata -S -l 26 --un $NAME.not_aligned_ensembl.fasta /nas3/gpratt/bowtie/ensembl_genes $NAME.not_alignednot_aligned_rRNA_all_ref.fasta $NAME.aligned.ensembl.sam

bowtie -m 15 --best -p 4 --strata -S -l 26 --un $NAME.not_aligned_hg19.fasta hg19 $NAME.not_aligned_ensembl.fasta $NAME.aligned.hg19.sam

samtools view -bS $NAME.aligned.ensembl.sam > $NAME.ensembl.bam
samtools sort $NAME.ensembl.bam $NAME.ensembl.sort
samtools index $NAME.emsembl.sort.bam

samtools view -bS $NAME.aligned.hg19.sam > $NAME.hg19.bam
samtools sort $NAME.hg19.bam $NAME.hg19.hg19.sort
samtools index $NAME.hg19.hg19.sort
