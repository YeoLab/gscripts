from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

try:
    name = sys.argv[1]+"_sort_bam"
except KeyError:
    name = "sort_bam"

cmd_list = []
for file in glob('*bam'):

	cmd_list.append('samtools sort -m 50000000000 {}  {}.sorted'.format(file, file))
sub = Submitter(queue_type='PBS', sh_file=name+'.sh', command_list=cmd_list, job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=16, queue='home', array=True, max_running=10)

