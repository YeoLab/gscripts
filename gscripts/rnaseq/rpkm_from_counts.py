from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

try:
    name = sys.argv[1]+"_rpkms_from_count"
except IndexError:
    name = "rpkms_from_count"

cmd_list = []
for file in glob('*count'):

	cmd_list.append('single_RPKM.py -i {} -o {}.rpkm'.format(file, file))

sub = Submitter(queue_type='PBS', sh_file=name+'.sh', command_list=cmd_list, job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', walltime='1:00:00', array=True, max_running=20)

