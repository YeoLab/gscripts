#!/usr/bin/env python

# Parse command line arguments
import argparse

# To submit jobs to the PBS queue
import gscripts.qtools as qtools

# Get filenames matching a pattern
from glob import glob

'''
Author: olga
Date created: 7/11/13 2:31 PM

The purpose of this program is to ...

Example run:
cd dir_of_sam_files
python sam_to_bam_and_sort.py
'''

# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='''Given a directory, convert all the sam files in it
             to bam, and sort the bam files
            ''',
            add_help=True, prefix_chars='-')
        self.parser.add_argument('--directory', '-d', action='store',
                                 type=str, default='',
                                 help='Where to look for sam files. Default '
                                      'is the current directory')
        self.parser.add_argument('--samtools-sort-args', action='store',
                                 type=str, default='', required=False,
                                 help='Arguments to pass to `samtools sort`, '
                                      'e.g. "-n" (needs to be in quotes) to '
                                      'sort by read name')
        self.parser.add_argument('job_name', required=False,
                                 type=str, action='store',
                                 default='repeat_align',
                                 help='Name of the job submitted to the cluster')

        self.parser.add_argument('--do-not-submit', required=False,
                                 action='store_true',
                                 help='Flag to not actually submit the job but '
                                      'just write the sh file (for testing)')
        self.parser.add_argument('-o', '--out-sh', required=False, type=str,
                                 action='store',
                                 help='Name of the sh file to submit to the job '
                                      'scheduler')
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


class SamToBamAndSort(object):
    def __init__(self, job_name, out_sh, directory, submit,
                 samtools_sort_args=''):
        # samtools_sort_args = cl.args['samtools_sort_args']
        # directory = cl.args['dir']

        # add final forward slash if it's not the current directory and it's
        # not the empty string. May cause bugs if your current directory is
        # the base directory '/', but I'm not too worried about that :)
        if directory != '' and not directory.endswith('/'):
            directory += '/'

        for sam in glob('%s*.sam' % directory):
            bam = sam.replace('.sam', '.bam')

            qsub_commands = []

            # Samtools view flags:
            # '-b': output is BAM
            # '-S': input is SAM
            qsub_commands.append('samtools view -bS %s > %s' % (sam, bam))

            sorted_prefix = bam.replace('.bam', '.sorted')
            qsub_commands.append('samtools sort %s %s %s' %
                                 (samtools_sort_args, bam, sorted_prefix))

            # submitter_prefix = 'sam2bam_sort_index_%s' % (bam)
            # submitter_sh = submitter_prefix + '.sh'
            sub = qtools.Submitter(queue_type='PBS',
                                   sh_filename=out_sh,
                                   commands=qsub_commands,
                                   job_name=job_name, nodes=1, ppn=16,
                                   queue='home')
            sub.write_sh(submit=submit)


# Function: main
def main():
    '''
    This function is invoked when the program is run from the command line,
    i.e. as:
        python program.py
    or as:
        ./program.py
    If the user has executable permissions on the user (set by chmod ug+x
    program.py or by chmod 775 program py. Just need the 4th bit set to true)
    '''

    try:
        cl = CommandLine()

        job_name = cl.args['job_name']
        out_sh = job_name + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']
        directory = cl.args['directory'].rstrip('/')

        SamToBamAndSort(job_name, out_sh, directory, submit,
                        cl.args['samtools_sort_args'])

    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()
