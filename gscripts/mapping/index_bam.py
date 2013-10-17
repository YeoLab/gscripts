from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

try:
    name = sys.argv[1]+"_index_bam"
except IndexError:
    name = "index_bam"

cmd_list = []
for file in glob('*sorted.bam'):

	cmd_list.append('samtools index {}'.format(file))
sub = Submitter(queue_type='PBS', sh_file=name+'.sh', command_list=cmd_list, job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', walltime='1:00:00', array=True, max_running=20)

