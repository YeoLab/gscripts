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
        self.parser.add_argument('--script-name', type=str, action='store',
                                 default='miso.sh',
                                 help='Assign a different name to this '
                                      'script, and the stderr/stdout produced'
                                      ' by the PBS job, too.')
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
            miso = (which('miso')[0])
            miso_scripts_dir = os.path.dirname(miso)
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
        num_processes = cl.args['num_processes']

        bams = []
        sample_ids = []
        notes = []

        with open(sample_info_file) as f:
            header = f.readline()
            # print 'header', header
            for line in f:
                # print 'line', line.rstrip().split('\t')
                sample_id, bam, note = line.rstrip().split('\t')
                sample_ids.append(sample_id)
                bams.append(bam)
                notes.append(note)

        # Command-line commands to submit to the cluster
        commands = []
        constitutive_exons_dir = '%s/%s_constitutive' % (base_annotation_dir, event_type)
        # print 'constitutive_exons_dir', constitutive_exons_dir
        # try:
        constitutive_exons_gff = glob('%s/*.gff' % constitutive_exons_dir)[0]

        commands.append('MISO=%s' % miso)
        commands.append('MISO_SCRIPTS_DIR=$(dirname $MISO)')
        commands.append('RUN_MISO_PY=$MISO_SCRIPTS_DIR/run_miso.py')
        commands.append('PAIRED_END_UTILS_PY=$MISO_SCRIPTS_DIR/pe_utils.py')
        commands.append('EVENT_TYPE=%s' % event_type)
        commands.append('EVENT_TYPE_INDEX=%s' % event_type_index)
        commands.append('CONSTITUTIVE_EXONS_GFF=%s' % constitutive_exons_gff)
        commands.append('BAMS_AND_IDS="%s"' %
                        ' '.join(','.join([bam, sample_id])
                                 for bam, sample_id in zip(bams, sample_ids) ))

        commands.append('EVENT_TYPE=%s\n' % event_type)
        commands.append('for i in $BAMS_AND_IDS ; do IFS=","')
        commands.append('    # Extract the data from the tuple')
        commands.append('    set $i')
        commands.append('    BAM=$1')
        commands.append('    ID=$2')
        commands.append('\n    echo')
        commands.append('    # Write down the sample ID and the time started '
                        'to know how long each one takes')
        commands.append('    echo "----- $ID -----"')
        commands.append('    date')
        commands.append('    DIR=$(dirname $BAM)')
        commands.append('    OUT_DIR=$DIR/miso/$EVENT_TYPE/$ID')
        commands.append("\n   # Create the output directory if it doesn't "\
                        "exist")
        commands.append('    if [ ! -d $OUT_DIR ] ; then')
        commands.append('        mkdir -p $OUT_DIR')
        commands.append('    fi')
        commands.append('\n    # Remove any *.bed files '\
                        'previously created by --prefilter')
        commands.append('    rm -rf $OUT_DIR/*.bed\n')
        commands.append("    # Create the insert size file if it doesn't exist"
                        "already")
        commands.append('    INSERT_SIZE_FILE=$BAM.insert_len')
        commands.append('    if [ ! -e $INSERT_SIZE_FILE ] ; then')
        commands.append('        INSERT_SIZE_COMMAND="python ' \
                        '$PAIRED_END_UTILS_PY $BAM $CONSTITUTIVE_EXONS_GFF '\
                        '$DIR"')
        commands.append('        date')
        commands.append('        echo Starting ... $INSERT_SIZE_COMMAND')
        commands.append('        $INSERT_SIZE_COMMAND')
        commands.append('    fi')
        commands.append('    # Extract the insert size mean and standard '
                        'deviation')
        commands.append("    INSERT_SIZE_MEAN=$(head -n 1 $INSERT_SIZE_FILE "
                            "| sed 's:#::' | cut -d'," \
                            "' -f1 | cut -d'=' -f2)")
        commands.append("    INSERT_SIZE_STDDEV=$(head -n 1 "
                            "$INSERT_SIZE_FILE "
                            "| sed 's:#::' | cut -d'," \
                            "' -f2 | cut -d'=' -f2)")
        commands.append('\n    # Assuming that the first read of the bam file'
                        ' is representative, such that all the reads in the '
                        '# file are exactly the same length, we can take the '
                        'first read from the bam file and measure its length,'
                        '# and use that for our algorithm')
        commands.append("    READ_LEN=$(samtools view $BAM | head -n 1 | "
                        "cut -f 10 | awk '{ print length }')")
        commands.append('\n    echo')
        commands.append('    date')
        commands.append('    MISO_COMMAND="python $MISO --run '
                        '$EVENT_TYPE_INDEX $BAM '
                        '--output_dir $OUT_DIR --read-len $READ_LEN '
                        '--paired-end '
                        '$INSERT_SIZE_MEAN $INSERT_SIZE_STDDEV -p %d '
                        '--no-filter-events"' % num_processes)
        commands.append('    echo Starting ...... $MISO_COMMAND')
        commands.append('    $MISO_COMMAND')
        commands.append('\n    # Now summarize the findings')
        commands.append('    '
                        'SUMMARIZE_COMMAND="python $RUN_MISO_PY '
                        '--summarize-samples $OUT_DIR $OUT_DIR"')
        commands.append('    echo')
        commands.append('    date')
        commands.append('    echo Starting .... $SUMMARIZE_COMMAND')
        commands.append('done')

        # Put the submitter script wherever the command was run from
        submit_sh = 'miso.sh'
        job_name = 'miso'
        sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                        command_list=commands, job_name=job_name)
        pbs_id = sub.write_sh(submit=True, nodes=1, ppn=16, queue='home-yeo')
        print pbs_id

        #TODO: write run comparisons script. skip duplicates
        # # Run comparisons
        # commands = []
        # commands.append('MISO=%s' % miso)
        # commands.append('MISO_SCRIPTS_DIR=$(dirname $MISO)')
        # commands.append('RUN_MISO_PY=$MISO_SCRIPTS_DIR/run_miso.py')
        # commands.append('PAIRED_END_UTILS_PY=$MISO_SCRIPTS_DIR/pe_utils.py')
        # commands.append('EVENT_TYPE=%s' % event_type)
        # commands.append('EVENT_TYPE_INDEX=%s' % event_type_index)
        # commands.append('CONSTITUTIVE_EXONS_GFF=%s' % constitutive_exons_gff)
        # commands.append('BAMS_AND_IDS="%s"' %
        #                 ' '.join(','.join([bam, sample_id])
        #                          for bam, sample_id in zip(bams, sample_ids) ))
        # commands.append('IDS="%s"' % ' '.join(sample_ids))
        # commands.append('EVENT_TYPE=%s\n' % event_type)
        #
        # commands.append('\nfor BAM,ID1 in $BAMS_AND_IDS ; do IFS=","')
        # commands.append('    DIR=$(dirname $BAM)')
        # commands.append('    mkdir -p $DIR/miso/$EVENT_TYPE/comparisons')
        # commands.append('    for ID2 in $IDS ; do')
        # commands.append('        if [ $ID1 != $ID2 ] ; then')
        # commands.append('            mkdir -p '
        #                 '$DIR/miso/$EVENT_TYPE/comparisons/$ID1\_vs_$ID2')
        # commands.append('            echo')


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