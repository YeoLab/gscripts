from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

try:
    name = sys.argv[1]+"_index_bam"
except IndexError:
    name = "index_bam"

command_list = []

for file in glob('*sorted.bam'):
    command_list.append('samtools index {0}'.format(file))


def submit_and_write(name, command_list):
    sub = Submitter(queue_type='PBS', sh_file=name + '.sh',
                    command_list=command_list,
                    job_name=name)

    sub.write_sh(submit=True, nodes=1, ppn=1, queue='home', array=True,
                 max_running=20, walltime='1:00:00')


## max number of jobs in an array on TSCC is 500
#if len(command_list) > 500:
#    command_list_list = [command_list[i:(i + 500)] for i in xrange(0,
#                                                                   len(
#                                                                       command_list),
#                                                                   500)]
#    for i, commands in enumerate(command_list_list):
#        submit_and_write('{}{}'.format(name, i + 1), commands)
#else:
submit_and_write(name, command_list)