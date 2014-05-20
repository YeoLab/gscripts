#!/usr/bin/env python

from glob import glob
from gscripts.qtools._Submitter import Submitter
import sys
import argparse


class CommandLine(object):
    """
    Check out the argparse documentation for more awesome things you can do!
    https://docs.python.org/2/library/argparse.html
    """

    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Index sorted bam files (required by some many '
                        'programs)')
        parser.add_argument('job_name', required=False,
                            type=str, action='store', default='index_bam',
                            help='Name of the job submitted to the cluster')
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


class SortBam(object):
    def __init__(self, job_name, out_sh, directory, submit):
        command_list = []
        for filename in glob('{}/*bam'.format(directory)):
            command_list.append(
                'samtools sort -@ 8 -m 50000000000 {0} {0}.sorted'
                .format(filename))


        # def submit_and_write(name, command_list):
        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=command_list,
                        job_name=job_name, nodes=1, ppn=8, queue='home',
                        array=True,
                        max_running=10, walltime='0:30:00')

        sub.write_sh(submit=submit)


if __name__ == '__main__':
    try:
        cl = CommandLine()

        job_name = cl.args['job_name']
        out_sh = job_name + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']
        directory = cl.args['directory'].rstrip('/')

        SortBam(job_name, out_sh, directory, submit)

    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)