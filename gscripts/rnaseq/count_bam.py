from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

orientation = sys.argv[1]
species = sys.argv[2]
try:
    name = sys.argv[3]+"_count_bam"
except IndexError:
    name = 'count_bam'

if species == 'hg19':
	annot_file = '/projects/ps-yeolab/genomes/hg19/gencode_v17/gencode.v17.annotation.exons.bed'
elif species == 'mm9':
	annot_file = '/projects/ps-yeolab/genomes/mm9/Mus_musculus.NCBIM37.64.fixed.exons.bed'
elif species == 'spikein':
    annot_file = '/projects/ps-yeolab/genomes/spikein/arraycontrol_spike_sequences.bed'

cmd_list = []
for file in glob('*sorted.bam'):

	cmd_list.append('count_tags.py --annotation_file {} -f {} -b {} -o {}.count'.format(annot_file,orientation, file, file))

sub = Submitter(queue_type='PBS', sh_file=name+'.sh', command_list=cmd_list, job_name=name)
sub.write_sh(submit=True, nodes=1, ppn=16, queue='home', walltime='1:00:00', array=True, max_running=20)

