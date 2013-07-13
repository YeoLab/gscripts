#!/usr/bin/env python

# Parse command line arguments
import argparse

# Submit jobs to the TSCC cluster
from qtools import Submitter

'''
Author: olga
Date created: 7/13/13 10:16 AM

The purpose of this program is to ...

Example run:
python /home/obotvinnik/gscripts/gscripts/rnaseq/submit_rnaseqc.py --base-dir \
/oasis/tscc/scratch/obotvinnik/gm12878/rna-seq \
--sample-file /home/obotvinnik/projects/alt_first_exon/gm12878_samples.txt \
--species hg19
'''

# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='''Perform RNA-Seq quality control pipeline RNA-SeQC
            on the samples described in the sample file.
            ''',
            add_help=True, prefix_chars='-')
        self.parser.add_argument('--rnaseqc-bin', action='store',
                                 type=str,
                                 help='The RNA-SeQC binary executable file to'
                                      ' use.',
                                 default='/home/yeo-lab/software/rnaseqc/RNA-SeQC_v1.1.7.jar')
        self.parser.add_argument('--base-dir', '-d', action='store',
                                 required=True, type=str,
                                 help='The base directory to use. The default'
                                      ' is to have the out directory to be '
                                      'this base directory + /rnaseqc, e.g.:'
                                      'if the base '
                                      'directory is ~/scratch/single_cell/,'
                                      'then the rseqc results will be in '
                                      '~/scratch/single_cell/rnaseqc')
        self.parser.add_argument('--sample-file', action='store',
                                 required=True, type=str,
                                 help='A tab-delimited description of the '
                                      'samples and their BAMs. The header is:'
                                      'Sample ID    Bam File    Notes')
        self.parser.add_argument('--species', '-s', action='store',
                                 required=True, type=str,
                                 help='Which species these reads were mapped '
                                      'to',
                                 default='hg19')
        self.parser.add_argument('--additional-arguments', action='store',
                                 required=False, type=str,
                                 default='',
                                 help='Any additional arguments to pass to '
                                      'RNA-SeQC, besides the required ones.')
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
        rnaseqc_bin = cl.args['rnaseqc_bin']
        species = cl.args['species']
        genome_fasta = '/projects/ps-yeolab/genomes/%s/chromosomes/%s.fa.fai' \
                       % (species, species)
        base_dir = cl.args['base_dir']

        # add trailing forward slash if it's not already there
        base_dir = base_dir if base_dir.endswith('/') else base_dir + '/'
        output_dir = base_dir + 'rseqc/'
        sample_file = cl.args['sample_file']
        additional_arguments = cl.args['additional_arguments']

        job_name = 'rnaseqc_%s' % base_dir
        submit_sh = '%srseqc.sh' % base_dir
        command = 'java -jar %s -o %s -r %s -s %s %s' % (rnaseqc_bin,
                                                         output_dir,
                                                         genome_fasta,
                                                         sample_file,
                                                         additional_arguments)
        commands = [command]

        submit_out = submit_sh + '.out'
        submit_err = submit_sh + '.err'

        sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                        command_list=commands, job_name=job_name)
        sub.add_resource('-o', submit_out)
        sub.add_resource('-e', submit_err)
        sub.write_sh(submit=True, nodes=1, ppn=16, queue='glean')
        pass
    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()