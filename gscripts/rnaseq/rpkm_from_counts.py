from glob import glob
from gscripts.qtools._Submitter import Submitter

cmd_list = []
for file in glob('*count'):

	cmd_list.append('single_RPKM.py -i {} -o {}.rpkm'.format(file, file))

sub = Submitter(queue_type='PBS', sh_file='rpkm_from_count.sh', command_list=cmd_list, job_name='rpkm_from_count')
sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', walltime='1:00:00', array=True, max_running=20)

