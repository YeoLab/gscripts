#!/usr/bin/env python
# coding=utf-8

# Parse command line arguments
import argparse
import os
import sys

# Submit jobs to the cluster
from ..qtools import Submitter


'''
Author: olga
Date created: 7/12/13 9:38 AM

The purpose of this program is to write submitter scripts to perform MISO
analysis on a large amount of files. This script assumes paired-end reads.

# **Note** for some reason,

Example runs:

submit_miso_pipeline.py --psi-only --paired-end --event-type SE \
--sample-info-file ~/projects/singlecell/singlecell/sample_info.txt \
--sh-scripts-dir ~/projects/singlecell/scripts/ \
--summary-output-dir-base ~/projects/singlecell/singlecell/analysis/ \
--psi-walltime '0:50:00' --annotation-index-strfmt \
'/home/yeo-lab/genomes/hg19/miso_annotations/%s_index' \
--sample-id-suffix min10reads \
--extra-miso-arguments ' --settings-filename \
~/MISO/misopy/settings/miso_settings_min_event_reads10.txt' \
--individual-jobs

submit_miso_pipeline.py --psi-only --paired-end --event-type SE --sample-info-file ~/projects/singlecell/singlecell/sample_info_round2and3.txt --sh-scripts-dir ~/projects/singlecell/scripts/ --summary-output-dir-base ~/projects/singlecell/singlecell/analysis/ --psi-walltime '1:50:00' --annotation-index-strfmt '/home/yeo-lab/genomes/hg19/miso_annotations/%s_index' --sample-id-suffix min10reads --extra-miso-arguments ' --settings-filename ~/MISO/misopy/settings/miso_settings_min_event_reads10.txt' --individual-jobs

submit_miso_pipeline.py --psi-only --single-end --event-type SE \
--sample-info-file ~/projects/singlecell/singlecell/sample_info.txt \
--sh-scripts-dir ~/projects/singlecell/scripts/ \
--summary-output-dir-base ~/projects/singlecell/singlecell/analysis/ \
--psi-walltime '0:50:00' \
--annotation-index-strfmt '/home/obotvinnik/genomes/hg19/miso_annotations/%s_index'


submit_miso_pipeline.py --run-all --single-bam $BAM $SAMPLE_ID\
--num-processes 16

TODO: deprecate the non-queue way of running
'''

# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='''Write a script to perform MISO analysis
            on individual samples. This script must be submitted to
            ''',
            add_help=True, prefix_chars='-')

        self.parser.add_argument('--sample-id', type=str,
                                 action='store',
                                 help='sample ID. required if using --bam',
                                 required=False)
        samples = self.parser.add_mutually_exclusive_group(required=True)
        samples.add_argument('--bam', type=str,
                             action='store', default='',
                             help='A single BAM file')
        samples.add_argument('--sample-info', type=str, action='store',
                             help='A tab-delimited file with Bam files as '
                                  'column 1 and ')

        self.parser.add_argument('--debug', action='store_true',
                                 default=False,
                                 help="Don't make any files, just print "
                                      "everything that would have been made")
        self.parser.add_argument('--genome', type=str, action='store',
                                 required=True, help='Which genome to use')
        self.parser.add_argument('--output-sh', type=str, required=True,
                                 action='store',
                                 help="The name of the .sh script created for one-touch action")
        self.parser.add_argument('--summary-output-dir-base', type=str,
                                 action='store', default='',
                                 help='The base directory for which to place '
                                      'the MISO summary outputs. By '
                                      'default, MISO outputs are of '
                                      'the format: (base_dir)/miso/('
                                      'event_type)/(sample_id). The default '
                                      'base dir is the directory of the .bam '
                                      'file, e.g. if the bam you provide is '
                                      'in ~/scratch/single_cell and your '
                                      'event type is "SE", then miso outputs '
                                      'for sample id "A1_02"'
                                      ' will be in the folder'
                                      '~/scratch/single_cell/miso/SE/A1_02/. '
                                      'However, the intermediate output would'
                                      ' still be in (bam_dir)/miso/('
                                      'event_type)/(sample_id) because MISO '
                                      'outputs a TON of intermediate files '
                                      'that nobody wants to deal with.'
                                      'Otherwise, if you provide a folder '
                                      'such as '
                                      '~/projects/single_cell/analysis, '
                                      'then the MISO summary output for the '
                                      'same '
                                      'sample would be in: '
                                      '~/projects/single_cell/analysis/miso/SE/A1_02')
        self.parser.add_argument('--sample-id-suffix', type=str,
                                 action='store', default='',
                                 help='Extra identification to add to these '
                                      'samples, e.g. if you are running with '
                                      'a settings file that specifies a '
                                      'minimum of 10 reads instead of 20, '
                                      'you could say "_min_event_reads10" as '
                                      'a suffix')
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
                                      'e.g. --extra-miso-arguments " '
                                      '--no-filter-events'
                                      ' --settings-filename '
                                      'miso_settings_min_event_reads5.txt". '
                                      'If this is not working for you, '
                                      'try adding a space between the first '
                                      'quote and the first dash of the miso '
                                      'argument. For some reason this helps..'
                                      '..')

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


class MisoPipeline(object):
    def __init__(self, bam, sample_id, output_sh,
                 genome,
                 read_type='single_end',
                 debug=False, num_processes=16,
                 submit_sh_suffix='',
                 sample_id_suffix='',
                 extra_miso_arguments=''):
        """
        Given a CommandLine object 'cl', save the arguments as
        attributes of this class. And initialize the cluster job IDs as None,
         so we can check if they're there or not in the future.
        """
        self.read_type = read_type
        #self.event_type = cl.args['event_type'].upper()
        # self.sample_info_file = cl.args['sample_info_file']

        self.debug = debug
        self.num_processes = num_processes

        self.bam = bam
        self.sample_id = sample_id
        self.sample_ids = [sample_id]
        self.bams = [bam]

        self.output_sh = output_sh

        self.extra_miso_arguments = extra_miso_arguments

        if submit_sh_suffix != '':
            self.submit_sh_suffix = '_' + submit_sh_suffix.lstrip('_')
        else:
            self.submit_sh_suffix = ''
        if sample_id_suffix != '':
            self.sample_id_suffix = '_' + sample_id_suffix.lstrip('_')
            self.sample_ids = [sample_id + self.sample_id_suffix
                               for sample_id in self.sample_ids]
        else:
            self.sample_id_suffix = ''

        #self.job_name_prefix = 'miso%s_%s' % (self.submit_sh_suffix,
        #                                      self.event_type)

        self.genome = genome


    def run_all_single_sample(self):

        bam = self.bam
        sample_id = self.sample_id

        bam_dir = os.path.dirname(os.path.abspath(bam))

        commands = []
        commands.append('#!/bin/bash')
        commands.append('# Finding all MISO splicing scores for sample: {}. '
                        'Yay!\n'
                        .format(sample_id))

        insert_len_arguments = ''

        event_types = ['SE', 'MXE', 'AFE', 'ALE', 'A3SS', 'A5SS',
                       'RI', 'TANDEMUTR']

        # Get the read length. Gonna keep this as bash because samtools
        # and less are very fast
        #commands.append(
        #    '\n# Assuming that the first read of the bam file is '
        #    'representative, such that all the reads in the '
        #    '\n# file are exactly the same length, we can take the first '
        #    'read from the bam file and measure its length, '
        #    '\n# and use that for our algorithm')
        commands.append(
            "READ_LEN=$(samtools view %s | head -n 1 | cut -f 10 | awk '{ "
            "print length }')" % (bam))

        for event_type in event_types:
            out_dir = '{}/miso/{}/{}'.format(os.path.dirname(os.path
                                                             .abspath(bam)),
                                             sample_id, event_type)
            psi_out = '{}/psi.out'.format(out_dir)
            psi_err = '{}/psi.err'.format(out_dir)

            commands.append('\n\n# calculate Psi scores for'
                            ' all {} events'.format(event_type))
            commands.append('mkdir -p {}'.format(out_dir))
            commands.append('python /projects/ps-yeolab/software/bin/miso \
 --run /projects/ps-yeolab/genomes/{0}/miso/{1}_index \
 {2} --output-dir {3} \
 --read-len $READ_LEN \
{4} \
 --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations'
                            '/miso_settings_min_event_reads10.txt \
 -p {5} \
 > {6} \
 2> {7}'.format(self.genome, event_type, bam, out_dir,
                insert_len_arguments, self.num_processes, psi_out,
                psi_err))

            commands.append("\n# Check that the psi calculation jobs didn't "
                            "fail.\n#'-z' "
                            "returns "
                            "true when a string is empty, so this is checking "
                            "that grepping these files for the words 'failed' "
                            "and 'shutdown' didn't find anything.")
            commands.append('iffailed=$(grep failed {})'.format(psi_out))
            commands.append('ifshutdown=$(grep shutdown {})'.format(psi_err))
            commands.append(
                "if [ ! -z \"$iffailed\" -o ! -z \"$ifshutdown\" ] ; "
                "then\n\
    #rm -rf {0}\n\
    echo \"MISO psi failed on event type: {1}\"\n\
    exit 1\n\
fi\n".format(out_dir, event_type))

            commands.append('# Summarize psi scores for all {} events'
                            .format(event_type))
            commands.append('python /home/yeo-lab/software/bin/run_miso.py '
                            '--summarize-samples {0} ' \
                            '{0} >{0}/summary.out 2>{0}/summary.err'.format(
                out_dir))
            commands.append("\n# Check that the summary jobs didn't fail")
            commands.append("# '-s' returns true if file size is nonzero, "
                            "and the error file should be empty.")
            commands.append("""if [ -s {0}/summary.err ] ; then
    #rm -rf {0}\n
    echo 'MISO psi failed on event type: {1}'
    exit 1
fi
""".format(out_dir, event_type))
        sh_file = self.output_sh
        with open(sh_file, 'w') as f:
            f.write('\n'.join(commands))
        sys.stdout.write('Wrote miso script for sample "{}": {}\n'.format(
            sample_id, sh_file))


# Function: main
def main():
    '''
    This function is invoked when the program is run from the command line,
    i.e. as:
        python submit_miso_pipeline.py
    or as:
        ./submit_miso_pipeline.py
    If the user has executable permissions on the user (set by chmod ug+x
    program.py or by chmod 775 program py. Just need the 4th bit set to true)
    '''
    cl = CommandLine()
    try:
        miso_pipeline = MisoPipeline(cl.args['bam'], cl.args['sample_id'],
                                     cl.args['output_sh'],
                                     cl.args['genome'],
                                     debug=cl.args['debug'],
                                     num_processes=cl.args['num_processes'],
                                     submit_sh_suffix=cl.args[
                                         'submit_sh_suffix'],
                                     sample_id_suffix=cl.args[
                                         'sample_id_suffix'],
                                     extra_miso_arguments=cl.args[
                                         'extra_miso_arguments'])

        # Read the arguments to see which piece of the MISO pipeline to run
        #if cl.args['run_all']:
        #    if cl.args['single_sample']:
        miso_pipeline.run_all_single_sample()
        #else:
        #    miso_pipeline.run_all()
        #elif cl.args['insert_len_only']:
        #    miso_pipeline.insert_len()
        #elif cl.args['psi_only']:
        #    miso_pipeline.psi()
        #elif cl.args['summary_only']:
        #    miso_pipeline.summary()
        #elif cl.args['psi_and_summary']:
        #    miso_pipeline.psi_and_summary()

    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()