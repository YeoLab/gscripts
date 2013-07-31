#!/usr/bin/env python
# coding=utf-8

# Parse command line arguments
import argparse

# Submit jobs to the cluster
from gscripts.qtools import Submitter

# File name manipulations
import os

# The program that does all the heavy lifting of creating the submitter scripts
from MisoPipeline import MisoPipeline


'''
Author: olga
Date created: 7/12/13 9:38 AM

The purpose of this program is to ...

Example run:
python ~/gscripts/gscripts/miso/submit_miso_pipeline.py --sample-info-file sample_info_individual_miso_failed.txt --event-type SE --submit-sh-suffix rerun_failed --psi-and-summary
'''

# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description=''' Given a number of digits "n" and number of
            iterations "N", calculate .....
            ''',
            add_help=True, prefix_chars='-')
        # self.parser.add_argument('--index-base-dir',
        #                          action='store',
        #                          type=str,
        #                          default='/home/obotvinnik/genomes/miso_annotations/hg19',
        #                          help='The base directory to use for '
        #                               'annotations. The annotation is assumed'
        #                               ' to be (index_base_dir)/('
        #                               'event_type)_indexed/')
        self.parser.add_argument('--event-type', '-e',
                                 action='store', type=str, required=True,
                                 help="Which event you'd like to index. One "
                                      "of:"+
                                      ('\n'
                                       '1. Skipped exons (SE)\n'
                                       '2. Alternative 3’/5’ splice sites ('
                                       'A3SS, A5SS)\n'
                                       '3. Mutually exclusive exons (MXE)\n'
                                       '4. Tandem 3’ UTRs (TandemUTR)\n'
                                       '5. Retained introns (RI)\n'
                                       '6. Alternative first exons (AFE)\n'
                                       '7. Alternative last exons (ALE)\n'
                                       '                                      '
                                      ) +
                                      "See http://genes.mit"
                                      ".edu/burgelab/miso/docs/#alternative-event-annotations for more information")
        self.parser.add_argument('--sample-info-file', required=True,
                                 type=str,
                                 action='store',
                                 help='A tab-delimited sample info file with '
                                      'the header:\n'
                                      'Sample ID\tBam File\t Notes')
        self.parser.add_argument('--miso-scripts-dir', type=str,
                                 action='store',
                                 help='Which directory to use as the prefix for '
                                      'miso scripts. Default is the directory'
                                      ' returned from the unix command line '
                                      'command "which miso".', required=False)
        self.parser.add_argument('--base-annotation-dir',
                                 type=str, action='store',
                                 help='Where the MISO annotations are housed.'
                                      ' The indexed version are assumed to be'
                                      ' [base_annotation_dir]/['
                                      'event_type]_index. For example, '
                                      'if the base annotation dir is '
                                      '/home/obotvinnik/genomes/miso_annotations/hg19 '
                                      'and the event type is AFE, '
                                      'then the annotations are assumed to be'
                                      ' in folder'
                                      '/home/obotvinnik/genomes/miso_annotations/hg19/AFE_index/',
                                 default='/home/obotvinnik/genomes/miso_annotations/hg19')
        # self.parser.add_argument('--read-len', '-l', type=int, action='store',
        #                          help='Read lengths. Assumed to be the same '
        #                               'for all samples', required=True)
        self.parser.add_argument('--num-processes', '-p', type=int,
                                 action='store', default=16,
                                 help='Number of subprocesses for MISO to run'
                                      '. If you are using a computing cluster'
                                      ' with several processors on a single '
                                      'node, use the number of processors '
                                      'you are requesting')
        self.parser.add_argument('--submit-sh-suffix', type=str,
                                 action='store',
                                 default='',
                                 help='Add a suffix to this '
                                      'script name, and the stderr/stdout '
                                      'produced'
                                      ' by the PBS job, too. The default is '
                                      'miso_[event_type].sh, for example if '
                                      'your event type is skipped exons (SE),'
                                      ' then the script is called miso_SE.sh '
                                      'If you add the argument '
                                      '"--submit-sh-suffix pooled" then the '
                                      'submit filename would be '
                                      '"miso_pooled_SE.sh"')
        self.parser.add_argument('--extra-miso-arguments', type=str,
                                 action='store',
                                 default='',
                                 help='Any additional MISO "compute psi" '
                                      'arguments you want'
                                      ' to supply to all the samples. The '
                                      'default is no additional arguments. '
                                      'Protect this argument with quotes so '
                                      'it does not get interpreted as an '
                                      'argument to the MISO pipeline script, '
                                      'e.g. --miso-arguments "--no-bam-filter'
                                      ' --settings-filename '
                                      'miso_settings_min_event_reads5.txt"')

        self.parser.add_argument('--sample-id-suffix', type=str,
                                 action='store', default='',
                                 help='Extra identification to add to these '
                                      'samples, e.g. if you are running with '
                                      'a settings file that specifies a '
                                      'minimum of 10 reads instead of 20, '
                                      'you could say "_min_event_reads10" as '
                                      'a suffix')

        self.parser.add_argument('--psi-walltime', type=str, action='store',
                                 default='24:00:00')
        self.parser.add_argument('--summary-walltime', type=str,
                                 action='store',
                                 default='24:00:00')

        # Which part of the pipeline do you want to run?
        pipeline_part = self.parser.add_mutually_exclusive_group(required=True)
        pipeline_part.add_argument('--insert-len-only',
                                 action='store_true', default=False,
                                 required=False,
                                 help='Only compute the insert lengths of the'
                                      ' bam files provided in the sample info'
                                      ' file. Need the event type for this '
                                      'because we will build the library of '
                                      'insert sizes from these events. A '
                                      'single job to the cluster will be '
                                      'submitted.')
        pipeline_part.add_argument('--psi-only',
                                   action='store_true', default=False,
                                   required=False,
                                   help='Only compute "psi" (percent-spliced-'
                                        'in) values for the provided samples.'
                                        ' Do not compute the insert lengths ('
                                        'assumes the insert length files are '
                                        'already there), '
                                        'and do not summarize. A single job '
                                        'to the cluster will be submitted.')
        pipeline_part.add_argument('--summary-only',
                                   action='store_true', default=False,
                                   help='Only compute the summary of all '
                                        '"psi" (percent-spliced-in) values '
                                        'for the provided samples, '
                                        'creating a tab-delimited file of '
                                        'every splicing event. Do not '
                                        'compute the insert lengths or the '
                                        'psi values themselves (assumes the '
                                        'psi values are already there). A '
                                        'single job to the cluster will be '
                                        'submitted.')
        pipeline_part.add_argument('--psi-and-summary',
                                   action='store_true', default=False,
                                   help='Compute the "psi" ('
                                        'percent-spliced-in) values for the '
                                        'provided samples, and summarize the '
                                        'output, which creates a '
                                        'tab-delimited file of every splicing'
                                        ' event. This is handy if you have '
                                        'already computed the insert lengths '
                                        'separately')
        pipeline_part.add_argument('--run-all', action='store_true',
                                   help='Compute the insert length mean and '
                                        'standard deviation, '
                                        'the "psi" (percent spliced-in) '
                                        'scores of all splicing events, '
                                        'and summarize the relevant events')

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
        miso_pipeline = MisoPipeline(cl)

        # Read the arguments to see which
        if cl.args['run_all']:
            miso_pipeline.run_all()
        elif cl.args['insert_len_only']:
            miso_pipeline.insert_len()
        elif cl.args['psi_only']:
            miso_pipeline.psi()
        elif cl.args['summary_only']:
            miso_pipeline.summary()
        elif cl.args['psi_and_summary']:
            miso_pipeline.psi_and_summary()




        '''
        ## Run MISO on a pair of paired-end sample (with insert length distribution with mean 250,
        ## standard deviation 15) using the mouse genome skipped exon annotations using the
        ## the cluster

        # Compute Psi values for control sample
        python run_events_analysis.py
        --compute-genes-psi mm9/pickled/SE data/control.bam --output-dir SE/control/
        --read-len 35 --paired-end 250 15 --use-cluster

        # Compute Psi values for knockdown sample
        python run_events_analysis.py
        --compute-genes-psi mm9/pickled/SE data/knockdown.bam --output-dir SE/knockdown/
         --read-len 35 --paired-end 250 15 --use-cluster


        ## Summarize the output (only run this once --compute-genes-psi finished!)
        ## This will create a "summary" directory in SE/control/ and in SE/knockdown/
        python run_miso.py --summarize-samples SE/control/ SE/control/
        python run_miso.py --summarize-samples SE/knockdown/ SE/knockdown/

        ## Detect differentially expressed isoforms between "control" and "knockdown"
        ## This will compute Bayes factors and delta Psi values between the samples
        ## and place the results in the directory SE/comparisons/control_vs_knockdown
        python run_miso.py --compare-samples SE/control/ SE/knockdown/ SE/comparisons/
        '''

    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()