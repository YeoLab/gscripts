__author__ = 'olga'

from glob import iglob, glob
from gscripts.qtools._Submitter import Submitter
import sys
import os
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
            description='Filter fastq files on quality score and others')
        parser.add_argument('job_name', required=False,
                            type=str, action='store', default='convert_sam',
                            help='Sample info file with bam files and sample '
                                 'IDs')
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


class FastqFilters(object):
    def __init__(self, job_name, out_sh, submit=False, directory='./'):
        try:
            os.mkdir('{}/filtered/'.format(directory.rstrip('/')))
        except OSError:
            pass

        commands = []
        for filename in iglob('{}/*.fastq.gz'.format(directory.rstrip('/'))):
            #TODO: the -l argument "20" should be a % of read length
            commands.append('echo {0}; zcat {0} | fastx_artifacts_filter | '
                            'fastq_quality_trimmer -l 20 -t 30 | '
                            'fastq_quality_filter -q 30 -p 90 -z '
                            '> filtered/{0}'.format(filename))

        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=commands,
                        job_name=job_name, nodes=1, ppn=2, queue='home',
                        array=True,
                        max_running=20, walltime='1:00:00')

        sub.write_sh(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()

    job_name = cl.args['job_name']
    out_sh = job_name + '.sh' if cl.args['out_sh'] is None \
        else cl.args['out_sh']
    submit = not cl.args['do_not_submit']
    directory = cl.args['directory'].rstrip('/')

    FastqFilters(job_name, out_sh, submit, directory)


