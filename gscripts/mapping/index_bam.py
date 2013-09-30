from glob import glob
from gscripts.qtools._Submitter import Submitter


cmd_list = []
for file in glob('*sorted.bam'):

	cmd_list.append('samtools index {}'.format(file))
sub = Submitter(queue_type='PBS', sh_file='index_bam.sh', command_list=cmd_list, job_name='index_bam')
sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', walltime='1:00:00', array=True, max_running=20)

