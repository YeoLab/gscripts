from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

orientation = sys.argv[1]
species = sys.argv[2]

if species == 'hg19':
	annot_file = '/home/gpratt/gencode.v17.exons.bed'
elif species == 'mm9':
	annot_file = '/home/gpratt/Mus_musculus.NCBIM37.64.fixed.exons.bed'


cmd_list = []
for file in glob('*sorted.bam'):

	cmd_list.append('count_tags.py --annotation_file {} -f {} -b {} -o {}.count'.format(annot_file,orientation, file, file))

sub = Submitter(queue_type='PBS', sh_file='count_bam.sh', command_list=cmd_list, job_name='count_bam')
sub.write_sh(submit=True, nodes=1, ppn=16, queue='home', walltime='1:00:00', array=True, max_running=20)

