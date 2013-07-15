#!/usr/bin/env python
# coding=utf-8

# Parse command line arguments
import argparse

# Submit jobs to the cluster
from gscripts.qtools import Submitter
#import gscripts
from gscripts import which

# File name manipulations
import os

# Exit the program cleanly
import sys

# Get groups of files with similar names
from glob import glob

'''
Author: olga
Date created: 7/12/13 9:38 AM

The purpose of this program is to ...

Example run:
python ~/gscripts/gscripts/miso/submit_miso_pipeline.py \
  --base-annotation-dir ~/genomes/miso_annotations/hg19/ \
  --event-type AFE \
  --sample-info-file ~/projects/alt_first_exon/gm12878_samples_read_groups.txt \
  --read-len 76
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
        self.parser.add_argument('--read-len', '-l', type=int, action='store',
                                 help='Read lengths. Assumed to be the same '
                                      'for all samples', required=True)
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
        event_type = cl.args['event_type']

        sample_info_file = cl.args['sample_info_file']

        try:
            miso_scripts_dir = os.path.dirname(which('miso')[0])
        except IndexError:
            # If there is an IndexError, that means that 'which' returned an
            # empty list, and thus there is no miso installed on the path.
            print >> sys.stderr, '"which miso" returned empty list, ' \
                                 'the program miso does not exist on your ' \
                                 'command line'
            sys.exit(1)

        run_events_analysis_py = '%s/run_events_analysis.py' % miso_scripts_dir
        paired_end_utils_py = '%s/pe_utils.py' % miso_scripts_dir
        base_annotation_dir = cl.args['base_annotation_dir']
        base_annotation_dir = base_annotation_dir if not base_annotation_dir\
            .endswith('/') else base_annotation_dir.rstrip('/')

        event_type_gff = glob('%s/%s*.gff' % (base_annotation_dir,
                                                    event_type))
        event_type_index = '%s/%s_index' % (base_annotation_dir, event_type)

        read_len = cl.args['read_len']


        bams = []
        sample_ids = []
        notes = []

        with open(sample_info_file) as f:
            header = f.readline()
            # print 'header', header
            for line in f:
                # print 'line', line
                sample_id, bam, note = line.rstrip().split('\t')
                sample_ids.append(sample_id)
                bams.append(bam)
                notes.append(note)

        # Command-line commands to submit to the cluster
        commands = []
        output_dirs = []

        for sample_id, bam, note in zip(sample_ids, bams, notes):
            # os.path.dirname returns the directory containing the bam file WITHOUT
            # the trailing forward slash: '/'
            # e.g.:
            # >>> os.path.dirname(
            # '/home/gpratt/projects/upf1/analysis/rna/318_UPF11_NoIndex_L004_R1.fq.polyATrim.adapterTrim.rmRep.sorted.bam')
            # '/home/gpratt/projects/upf1/analysis/rna'
            base_bam_dir = os.path.dirname(bam)
            output_dir = '%s/miso/%s' % (base_bam_dir, event_type)
            output_dirs.append(output_dir)
            try:
                os.makedirs(output_dir)
            except:
                # This is just to silence the error of directory creation
                # if it's
                pass

            # Get the insert size from the *.insert_len file created by
            # pe_utils.py
            # First line of '.insert_len' file:
            # #mean=158.4,sdev=13.7,dispersion=1.1,num_pairs=20091
            # This file is assumed to be the '.bam' file + '.insert_len'
            insert_size_file = bam + '.insert_len'

            # first compute insert sizes
            # python ~/obot_virtualenv/bin/pe_utils.py --compute-insert-len $(
            # echo $BAMS | tr ' ' ,) ~/genomes/miso_annotations/hg19/AFE_constitutive/AFE.hg19.min_20.const_exons.gff --no-bam-filter --output-dir insert_sizes

            try:
                with open(insert_size_file) as f:
                    pass

            except IOError:
                # test that the corresponding constitutive exon file exists
                constitutive_exons_dir = '%s%s_constitutive' % (
                    base_annotation_dir, event_type)
                try:
                    constitutive_exons_gff = glob('%s/*.gff' %
                                              constitutive_exons_dir)[0]
                    with open(constitutive_exons_gff) as f:
                        pass
                except IndexError:
                    # Make the constitutive exons gff file for finding
                    exon_utils = '%s/exon_utils.py' % miso_scripts_dir
                    event_type_constitutive_dir = '%s/%s_constitutive/' \
                                                  % (base_annotation_dir,
                                                     event_type)
                    exon_utils_command = 'python %s --get-const-exons %s ' \
                                 '--output-dir %s' \
                              % (exon_utils, event_type_gff,
                                 event_type_constitutive_dir)
                    commands.append(exon_utils_command)

                    # Make sure the exon_utils.py commmand of finding
                    # constitutive exons finished before continuing on to
                    # find the insert length mean and standard deviation
                    commands.append('sleep 500')

                constitutive_exons_gff = glob('%s/*.gff' %
                                              constitutive_exons_dir)[0]
                paired_end_utils_command = 'python %s --compute-insert-len ' \
                                           '%s %s --no-bam-filter ' \
                                           '--output-dir %s' \
                                           % (paired_end_utils_py, bam,
                                              constitutive_exons_gff,
                                              base_bam_dir)
                commands.append(paired_end_utils_command)

                # Make sure the insert size calculation finishes before the
                # run_events_analysis.py command
                commands.append('sleep 500')

            with open(insert_size_file) as f:
                # remove the starting comment mark '#' and split on commas
                line = f.readline().lstrip('#').split(',')

                # Keep the insert size mean and standard dev as a string so
                # we don't have to deal with conversion
                insert_size_mean = line[0].split('=')[1]
                insert_size_stddev = line[1].split('=')[1]

            run_events_analysis_command = 'python %s' \
                      ' --compute-genes-psi %s %s --output-dir %s' \
                      ' --read-len %d --paired-end %s %s' \
                      % (run_events_analysis_py, event_type_index, bam,
                         output_dir, read_len, insert_size_mean,
                         insert_size_stddev)
            commands.append(run_events_analysis_command)

        for output_dir in output_dirs:
            run_miso_py = '%s/run_miso.py' % miso_scripts_dir
            summarize_command = 'python %s --summarize-samples %s %s' \
                                % (run_miso_py, output_dir, output_dir)
            commands.append(summarize_command)

        # Put the submitter script wherever the command was run from
        submit_sh = 'submit_miso.sh'
        job_name = 'miso_pipeline'
        sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                        command_list=commands, job_name=job_name)
        print sub.write_sh(submit=True, nodes=16, ppn=2, queue='glean')
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