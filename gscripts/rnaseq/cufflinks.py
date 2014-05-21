__author__ = 'olga'

from glob import iglob, glob
from gscripts.qtools._Submitter import Submitter
import sys
import os
from gscripts.qtools import Submitter
import argparse
import os


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Calculate expression using cufflinks')
        parser.add_argument('orientation', required=True, action='store',
                            type=str,
                            help='The orientation you want to count. One of '
                                 '"both" (unstranded), "flip", (stranded data)')
        annotation = parser.add_mutually_exclusive_group(required=True)
        annotation.add_argument('species', type=str, action='store',
                                help='The species you want to do count_bam on'
                                     '. This will choose a reasonable default '
                                     '.bed file to count on')
        annotation.add_agrument('--gtf', type=str, action='store',
                                help='The .gtf file you want to count over')
        parser.add_argument('name', type=str, action='store', required=False,
                            help='Prefix for the name of the job to submit')
        parser.add_argument('--out-sh', action='store', type=str,
                            required=False,
                            help='The sh file written and submitted to the '
                                 'cluster')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true', default=False,
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')
        parser.add_argument('--queue-type', required=False, type=str,
                            action='store', default='PBS',
                            help='Type of the queue to submit to. For testing '
                                 'purposes on non-server devices, e.g. laptops')
        parser.add_argument('-d', '--directory', required=False,
                            action='store', default='./',
                            help='Directory where the bam files are. Default '
                                 'is the current directory.')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))

    def do_usage_and_die(self, str):
        '''
        If a critical error is encountered, where it is suspected that the
        program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution
        of the program.
        '''
        import sys

        print >> sys.stderr, str
        self.parser.print_usage()
        return 2


# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


class Cufflinks(object):
    def __init__(self, gtf, job_name, out_sh, submit, directory):
        commands = []
        for bam in iglob('{}/*.sorted.bam'.format(directory.rstrip('/'))):
            commands.append('cufflinks --GTF {0} --GTF-guide  '
                            '--multi-read-correct --num-threads 8 {1}'.format(
                gtf, bam
            ))

        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=commands,
                        job_name=job_name, nodes=1, ppn=8, walltime='0:30:00',
                        array=True,
                        max_running=20
        )
        sub.write_sh(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()
    try:
        job_name = '_'.join([cl.args['name'], 'count_bam'])

        out_sh = name = job_name + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']
        directory = cl.args['directory']

        species = cl.args['species']

        if species == 'hg19':
            gtf = '/projects/ps-yeolab/genomes/hg19/gencode_v17/gencode.v17.annotation.gtf'
        elif species is None:
            if cl.args['gtf'] is not None:
                gtf = cl.args['gtf']
            else:
                raise ValueError('No known species provided, or alternative '
                                 'gtf provided')
        Cufflinks(gtf, job_name, out_sh, submit, directory)


    except Usage, err:
        cl.do_usage_and_die()

