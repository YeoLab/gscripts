#!/usr/bin/env python

# Parse command line arguments
import argparse

'''
Author: olga
Date created: 7/12/13 9:38 AM

The purpose of this program is to ...

Example run:
python submit_miso_pipeline.py -N 10
'''

# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description=''' Given a number of digits "n" and number of
            iterations "N", calculate .....
            ''',
            add_help=True, prefix_chars='-')
        self.parser.add_argument('--index-base-dir',
                                 action='store',
                                 type=str,
                                 default='/home/obotvinnik/genomes/miso_annotations/hg19',
                                 help='The base directory to use for '
                                      'annotations. The annotation is assumed'
                                      ' to be (index_base_dir)/('
                                      'event_type)_indexed/')
        self.parser.add_argument('--event-type', '-e',
                                 action='store', type=str, required=True,
                                 help="Which event you'd like to index. One "
                                      "of: See http://genes.mit"
                                      ".edu/burgelab/miso/docs/#alternative-event-annotations for more information")
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
        N = cl.args['N']

## Run MISO on a pair of paired-end sample (with insert length distribution with mean 250,
## standard deviation 15) using the mouse genome skipped exon annotations using the
## the cluster

# Compute Psi values for control sample
python run_events_analysis.py --compute-genes-psi mm9/pickled/SE data/control.bam --output-dir SE/control/ --read-len 35 --paired-end 250 15 --use-cluster

# Compute Psi values for knockdown sample
python run_events_analysis.py --compute-genes-psi mm9/pickled/SE data/knockdown.bam --output-dir SE/knockdown/ --read-len 35 --paired-end 250 15 --use-cluster


## Summarize the output (only run this once --compute-genes-psi finished!)
## This will create a "summary" directory in SE/control/ and in SE/knockdown/
python run_miso.py --summarize-samples SE/control/ SE/control/
python run_miso.py --summarize-samples SE/knockdown/ SE/knockdown/

## Detect differentially expressed isoforms between "control" and "knockdown"
## This will compute Bayes factors and delta Psi values between the samples
## and place the results in the directory SE/comparisons/control_vs_knockdown
python run_miso.py --compare-samples SE/control/ SE/knockdown/ SE/comparisons/

        pass
    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()