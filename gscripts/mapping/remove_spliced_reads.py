__author__ = 'olga'
from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

name = "remove_spliced"

try:
    name = sys.argv[1] + name
except IndexError:
    pass

cmd_list = []
for file in glob('*bam'):
    cmd_list.append("samtools view -h -F 4 {0} | awk '$6 !~ /N/ || $1 ~ /@/' "
                    "| "
                    "samtools view -bS - > {0}.unspliced.bam".format(file))

sub = Submitter(queue_type='PBS', sh_file=name + '.sh', command_list=cmd_list,
                job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=16, queue='home', array=True,
             walltime='0:30:00',
             max_running=10)
