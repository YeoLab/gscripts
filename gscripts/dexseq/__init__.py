__author__ = 'olga'

'''
How to "prepare" a DEXSeq annotation:
Basically:
python dexseq_prepare_annotation.py <in.gtf> <out.gff>

Specifically:
python /home/yeo-lab/software/dexseq_scripts/dexseq_prepare_annotation.py \
    /projects/ps-yeolab/genomes/hg19/gencode_v17/gencode.v17.annotation.gtf \
    ~/genomes/hg19/gencode_v17/gencode.v17.annotation.dexseq.gff


How to run DEXSeq:
python dexseq_count.py [options] <flattened_gff_file> <sam_file> <output_file>

From:
(obot_virtualenv)[obotvinnik@tscc-login2 single_cell]$ python $(which dexseq_count.py)
Usage: python dexseq_count.py [options] <flattened_gff_file> <sam_file> <output_file>

This script counts how many reads in <sam_file> fall onto each exonic part
given in <flattened_gff_file> and outputs a list of counts in <output_file>,
for further analysis with the DEXSeq Bioconductor package. (Notes: The
<flattened_gff_file> should be produced with the script
dexseq_prepare_annotation.py). <sam_file> may be '-' to indicate standard
input.

Options:
  -h, --help            show this help message and exit
  -p PAIRED, --paired=PAIRED
                        'yes' or 'no'. Indicates whether the data is paired-
                        end (default: no)
  -s STRANDED, --stranded=STRANDED
                        'yes', 'no', or 'reverse'. Indicates whether the data
                        is from a strand-specific assay (default: yes ). Be
                        sure to switch to 'no' if you use a non strand-
                        specific RNA-Seq library preparation protocol.
                        'reverse' inverts strands and is neede for certain
                        protocols, e.g. paired-end with circularization.
  -a MINAQUAL, --minaqual=MINAQUAL
                        skip all reads with alignment quality lower than the
                        given minimum value (default: 10)

Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology
Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General
Public License v3. Part of the 'DEXSeq' package.
'''
