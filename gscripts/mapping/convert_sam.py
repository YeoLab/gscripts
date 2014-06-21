#!/usr/bin/env python

import argparse
from glob import glob
import sys
from gscripts.qtools._Submitter import Submitter


class CommandLine(object):
    """
    Check out the argparse documentation for more awesome things you can do!
    https://docs.python.org/2/library/argparse.html
    """

    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Convert same files to bam files (binary format, '
                        'better compression). Also most programs take .bam '
                        'files as input')
        parser.add_argument('--job_name', required=False,
                            type=str, action='store', default='convert_sam',
                            help='Name of submitted job to scheduler')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true',
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')
        parser.add_argument('-o', '--out-sh', required=False, type=str,
                            action='store',
                            help='Name of the sh file to submit to the job '
                                 'scheduler')
        parser.add_argument('--queue-type', required=False, type=str,
                            action='store', default='PBS',
                            help='Type of the queue to submit to. For testing '
                                 'purposes on non-server machines (e.g. '
                                 'laptops)')
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


class ConvertSam(object):
    def __init__(self, job_name, out_sh=None, queue_type='PBS',
                 directory='./', submit=True):
        cmd_list = []
        for file in glob('{}/*sam'.format(directory)):
            cmd_list.append(
                'samtools view -bS -q 10 {} > {}.bam'.format(file, file))

        sub = Submitter(queue_type=queue_type,
                        sh_filename=out_sh,
                        commands=cmd_list,
                        job_name=job_name, nodes=1, ppn=1, queue='home',
                        walltime='1:00:00',
                        array=True, max_running=20)
        sub.job(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()

    job_name = cl.args['job_name']
    out_sh = job_name + '.sh' if cl.args['out_sh'] is None \
        else cl.args['out_sh']
    submit = not cl.args['do_not_submit']
    directory = cl.args['directory'].rstrip('/')

    ConvertSam(job_name=job_name, out_sh=out_sh,
               queue_type=cl.args['queue_type'], directory=directory,
               submit=submit)



