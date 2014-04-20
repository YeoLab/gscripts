#!/usr/bin/env python

from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys
import os

species = sys.argv[1]
try:
    jobname = sys.argv[2] + "_map_to_" + species
except IndexError:
    jobname = "map_to_" + species

cmd_list = []

# for file in glob('*R1*fastq'):
#     pair = file.replace('R1', 'R2')
#     name = file.replace('_R1', '')
#     cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
# --runMode alignReads \
# --runThreadN 16 \
# --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
# --genomeLoad LoadAndRemove \
# --readFilesIn {} {} \
# --outFileNamePrefix {}. \
# --outSAMunmapped Within \
# --outFilterMultimapNmax 1'.format(species, file, pair, name))
#
# for file in glob('*R1*norep'):
#     pair = file.replace('R1', 'R2')
#     name = file.replace('_R1', '')
#     cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
# --runMode alignReads \
# --runThreadN 16 \
# --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
# --genomeLoad LoadAndRemove \
# --readFilesIn {} {} \
# --outFileNamePrefix {}. \
# --outSAMunmapped Within \
# --outFilterMultimapNmax 1'.format(species, file, pair, name))

pwd = os.path.abspath(os.path.curdir)

for read1 in glob('*R1*gz'):
    read2 = read1.replace('R1', 'R2')
    name = '_'.join(read1.split('_'[:2]))
    cmd_list.append('STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir /projects/ps-yeolab/genomes/{0}/star_sjdb/ \
--genomeLoad LoadAndRemove \
--readFilesCommand zcat \
--readFilesIn {1}/{2} {1}/{3} \
--outFileNamePrefix {1}/{4}. \
--outSAMunmapped Within \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 5 \
--clip5pNbases 10 \
--clip3pNbases 10 \
--outFilterMultimapNmax 5'.format(species, pwd, read1, read2, name))

sub = Submitter(queue_type='PBS', sh_file=jobname + '.sh',
                command_list=cmd_list,
                job_name=jobname)
sub.write_sh(submit=True, nodes=1, ppn=8, walltime='0:30:00', use_array=True,
             array=True,
             max_running=20)
