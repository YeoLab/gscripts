#!/usr/bin/env python

from glob import iglob, glob
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

sample_ids = set([])

for read1 in iglob('*R1*gz'):
    sample_id = '_'.join(read1.split('_')[:2])
    if sample_id in sample_ids:
        continue

    read1 = ','.join(glob('{}*R1*gz'.format(sample_id)))
    read2 = read1.replace('R1', 'R2')
    sample_ids.add(sample_id)

    # print sample_id
    cmd_list.append('STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir /projects/ps-yeolab/genomes/{0}/star_sjdb/ \
--genomeLoad LoadAndRemove \
--readFilesCommand zcat \
--readFilesIn {2} {3} \
--outFileNamePrefix aligned/{4}. \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 5 \
--clip5pNbases 10 \
--clip3pNbases 10 \
--outFilterScoreMin 10 \
--outSAMattributes All \
--outFilterMultimapNmax 5'.format(species, pwd, read1, read2, sample_id))

sub = Submitter(queue_type='PBS', sh_file=jobname + '.sh',
                command_list=cmd_list,
                job_name=jobname)
sub.write_sh(submit=True, nodes=1, ppn=8, walltime='0:20:00', use_array=True,
             array=True,
             max_running=20)
