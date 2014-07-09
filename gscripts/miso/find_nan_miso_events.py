#!/usr/bin/env python

__author__ = 'olga'

from gscripts.qtools import Submitter
import argparse
from glob import iglob
import os


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Search miso outputs for log_score of nan and save '
                        'them to a file for that sample')
        parser.add_argument('--name', required=False,
                            type=str, action='store',
                            default='find_nan_miso_events',
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
                            help='Base directory, which has a "miso" directory there. This assumes '
                                 'the following directory structure:'
                                 '\n<directory>/miso/<sample_id>/<event_type>'
                                 '\nWhere "<directory>" is the location '
                                 'specified through this '
                                 'variable. If you ran your MISO samples using the Yeo Lab '
                                 'pipeline, you are fine. Default '
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


class FindNanMisoEvents(object):
    def __init__(self, job_name, out_sh, directory='./', queue_type='PBS',
                 submit=False):
        """
        Given a base folder which has a directory called "miso" where all the
        miso output is, search for bad events in the subfolders and then
        write them to a "nan_events.txt" file for that sample and event type.

        Parameters
        ----------
        job_name : str
            Name of the array job to be submitted
        out_sh : str
            Filename to write all the submitter commands to
        directory : str
            Base directory, which has a "miso" directory there. This assumes
            the following directory structure:
            <directory>/miso/<sample_id>/<event_type>
            Where "<directory>" is the location specified through this
            variable. If you ran your MISO samples using the Yeo Lab
            pipeline, you're fine.

        """
        commands = []

        glob_command = '{}/miso/*/*'.format(directory.rstrip('/'))
        for folder in iglob(glob_command):
            command = 'cd %s ; grep -m 1 nan chr*/*.miso > nan_events.txt' % \
                      folder
            commands.append(command)

        sub = Submitter(queue_type=queue_type, job_name=job_name,
                        sh_filename=out_sh,
                        commands=commands,
                        nodes=1, ppn=1, queue='home',
                        array=True,
                        max_running=20,
        )
        sub.write_sh(submit=submit)


if __name__ == '__main__':
    try:
        cl = CommandLine()

        job_name = cl.args['name']
        out_sh = name = job_name + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']
        directory = cl.args['directory'].rstrip('/')
        queue_type = cl.args['queue_type']

        FindNanMisoEvents(job_name, out_sh, directory, queue_type, submit)

    except Usage, err:
        cl.do_usage_and_die()