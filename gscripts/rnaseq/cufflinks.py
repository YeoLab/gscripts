__author__ = 'olga'

from glob import iglob, glob
from gscripts.qtools._Submitter import Submitter
import sys
import os

species = sys.argv[1]
try:
    jobname = sys.argv[2] + "_cufflinks_" + species
except IndexError:
    jobname = "cufflinks_" + species

commands = []

if species == 'hg19':
    gtf = '/projects/ps-yeolab/genomes/hg19/gencode_v17/gencode.v17.annotation.gtf'

for bam in iglob('*.sorted.bam'):
    commands.append('cufflinks --GTF {0} --GTF-guide  '
                    '--multi-read-correct --num-threads 8 {1}'.format(
        gtf, bam
    ))

sub = Submitter(queue_type='PBS', sh_file=jobname + '.sh',
                command_list=commands,
                job_name=jobname)
sub.write_sh(submit=True, nodes=1, ppn=8, walltime='0:30:00', use_array=True,
             array=True,
             max_running=20)