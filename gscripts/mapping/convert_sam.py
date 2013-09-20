from glob import glob
from gscripts.qtools._Submitter import Submitter

cmd_list = []
for file in glob('*sam'):

	cmd_list.append('samtools view -bS {} > {}.bam'.format(file, file))
sub = Submitter(queue_type='PBS', sh_file='convert_sam.sh', command_list=cmd_list, job_name='convert_sam')
sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', walltime='1:00:00', array=True, max_running=20)

