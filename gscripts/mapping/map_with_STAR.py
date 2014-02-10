#!/usr/bin/env python

from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

species = sys.argv[1]
try:
    name = sys.argv[2]+"_map_to_"+species
except IndexError:
    name = 'map_to_'+species
cmd_list = []
for file in glob('*fastq'):

	cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir /projects/ps-yeolab/genomes/{}/star/ \
--genomeLoad LoadAndRemove \
--readFilesIn {} --outFileNamePrefix {}. \
--outSAMunmapped Within \
--outFilterMultimapNmax 10'.format(species, file, file))

for file in glob('*txt'):

	cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir /projects/ps-yeolab/genomes/{}/star/ \
--genomeLoad LoadAndRemove \
--readFilesIn {} --outFileNamePrefix {}. \
--outSAMunmapped Within \
--outFilterMultimapNmax 10'.format(species, file, file))

for file in glob('*norep'):

	cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir /projects/ps-yeolab/genomes/{}/star/ \
--genomeLoad LoadAndRemove \
--readFilesIn {} --outFileNamePrefix {}. \
--outSAMunmapped Within \
--outFilterMultimapNmax 10'.format(species, file, file))

for file in glob('*gz'):

	cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir /projects/ps-yeolab/genomes/{}/star/ \
--genomeLoad LoadAndRemove \
--readFilesCommand zcat \
--readFilesIn {} --outFileNamePrefix {}. \
--outSAMunmapped Within \
--outFilterMultimapNmax 10'.format(species, file, file))

sub = Submitter(queue_type='PBS', sh_file=name+'.sh', command_list=cmd_list, job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=16, walltime='2:00:00', array=True, max_running=20)
