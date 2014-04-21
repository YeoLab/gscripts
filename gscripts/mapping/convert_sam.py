#!/usr/bin/env python


from glob import glob
import sys
from gscripts.qtools._Submitter import Submitter

try:
    name = sys.argv[1] + '_convert_sam'
except IndexError:
    name = 'convert_sam'

cmd_list = []
for file in glob('*sam'):
    cmd_list.append('samtools view -bS -q 10 {} > {}.bam'.format(file, file))

sub = Submitter(queue_type='PBS', sh_file=name + '.sh', command_list=cmd_list,
                job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', walltime='0:30:00',
             array=True, max_running=20)

