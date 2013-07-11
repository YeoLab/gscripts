__author__ = 'olga'

from glob import glob
from qtools import Submitter


hg19 = '/projects/ps-yeolab/genomes/hg19/'

genome_base_dir = '/projects/ps-yeolab/genomes/'

STAR = '/home/yeo-lab/software/bin/STAR'

cmd_list = ['%s ']
sub = Submitter(queue_type='PBS', sh_file='index_' + file + '.sh',
                command_list=cmd_list, job_name='index_' + file)
sub.write_sh(submit=True, nodes=1, ppn=1, queue='glean')