__author__ = 'olga'

from glob import iglob, glob
from gscripts.qtools._Submitter import Submitter
import sys
import os

try:
    os.mkdir('filtered/')
except OSError:
    pass

commands = []
for filename in iglob('*.fastq.gz'):
    #TODO: the -l argument "20" should be a % of read length
    commands.append('zcat {0} | fastx_artifacts_filter | '
                    'fastq_quality_trimmer -l 20 -t 30 | '
                    'fastq_quality_filter -q 30 -p 90 -z '
                    '> filtered/{0}'.format(filename))


def submit_and_write(name, command_list):
    sub = Submitter(queue_type='PBS', sh_file=name + '.sh',
                    command_list=command_list,
                    job_name=name)

    sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', array=True,
                 max_running=20, walltime='0:30:00')


submit_and_write('fastx_filter', commands)