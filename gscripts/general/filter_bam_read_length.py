#!/usr/bin/env python

from __future__ import print_function

# Parse command line arguments
import argparse

import pysam

'''
Author: Olga Botvinnik, based of off Emily Wheeler's code
Date created: 9/14/16 10:51 PM

The purpose of this program is to filter bam files for reads of a certain length

Example run:
filter_bam_read_length.py --length 100 aligned.bam
'''


# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='''Filter a bam file for only reads of a certain length''',
            add_help=True, prefix_chars='-')
        self.parser.add_argument('bam', action='store',
                                 type=str, help='Alignment file to filter')
        self.parser.add_argument('--length', '-l', action='store',
                                 type=int, default=100, required=False,
                                 help='Length of reads to filter for '
                                      '(default=100)')
        output_parser = self.parser.add_mutually_exclusive_group(
            required=False)
        output_parser.add_argument('--suffix', '-s', action='store',
                                   type=str, default='.100bp.bam',
                                   required=False,
                                   help='What to replace the ".bam" in the'
                                        ' file extension with. Not compatible '
                                        'with --output. '
                                        '(default=".100bp.bam")')
        output_parser.add_argument('--output', '-o', action='store', type=str,
                                   default=None, required=False,
                                   help='The full name of the new bam file to '
                                        'create. Default is to create a new '
                                        'bam file with the suffix .100bp.bam '
                                        'in the same folder (aka use the '
                                        '--suffix command). Cannot be used '
                                        'with the --suffix command.')
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


class FilterBamOnReadLength(object):

    def __init__(self, bam, length, suffix, output=None):
        if output is None:
            outfilename = bam.replace('.bam', suffix)
        else:
            outfilename = output

        infile = pysam.AlignmentFile(bam, "rb")
        outfile = pysam.AlignmentFile(outfilename, 'wb', template=infile)

        print('Getting only reads in {bam} that are exactly '
                         '{length}bp ...\n'.format(bam=bam, length=length))
        for read in infile.fetch():
            if read.infer_query_length() == length:
                outfile.write(read)
        print('\tDone.')

        outfile.close()
        infile.close()

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
    cl = CommandLine()
    try:
        FilterBamOnReadLength(**cl.args)
    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()