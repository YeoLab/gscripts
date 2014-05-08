#!/usr/bin/env python

from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys

try:
    name = sys.argv[1] + "_sort_bam"
except IndexError:
    name = "sort_bam"

command_list = []
for file in glob('*bam'):
    command_list.append('samtools sort -@ 8 -m 50000000000 {0} {0}.sorted'
                        .format(
        file))


def submit_and_write(name, command_list):
    sub = Submitter(queue_type='PBS', sh_file=name + '.sh',
                    command_list=command_list,
                    job_name=name)

    sub.write_sh(submit=True, nodes=1, ppn=8, queue='home', array=True,
                 max_running=10, walltime='0:30:00')


# max number of jobs in an array on TSCC is 500
#if len(command_list) > 500:
#    command_list_list = [command_list[i:(i + 500)] for i in xrange(0,
#                                                                   len(
#                                                                       command_list),
#                                                                   500)]
#    for i, commands in enumerate(command_list_list):
#        submit_and_write('{}{}'.format(name, i + 1), commands)
#else:
submit_and_write(name, command_list)
