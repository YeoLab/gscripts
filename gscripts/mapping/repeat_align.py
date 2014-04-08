#!/usr/bin/env python

from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

cmd_list = []
for file in glob('*fastq'):
    cmd_list.append('bowtie \
-c \
-S \
-q \
-p 16 \
-e 100 \
-l 20 \
--un {}.norep \
all_ref \
{} \
| grep -v \"@\" \
|  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl \
> {}.repeat_counts'.format(file, file, file))

for file in glob('*gz'):
    cmd_list.append('gunzip -c {} \
|bowtie \
-c \
-S \
-q \
-p 16 \
-e 100 \
-l 20 \
--un {}.norep \
all_ref \
- \
| grep -v \"@\" \
|  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl \
> {}.repeat_counts'.format(file, file, file))

sub = Submitter(queue_type='PBS', sh_file='repeat_align.sh',
                command_list=cmd_list, job_name='repeat_align')
sub.write_sh(submit=True, nodes=1, ppn=16, walltime='2:30:00', array=True,
             max_running=20)
