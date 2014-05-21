__author__ = 'olga'

from gscripts.qtools import Submitter
import os
from glob import glob
import sys
from gscripts import which
import pysam
from gscripts.general import read_sample_info_file


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
        # self.sample_ids = [sample_id]
        # self.bams = [bam]

        self.output_sh = output_sh

        self.extra_miso_arguments = extra_miso_arguments

        # Initialize the IDs we're going to use
        self.insert_len_job_id = None
        #self.insert_len_job_id = dict((sample_id, None) for sample_id in
        #                              self.sample_ids)
        #self.psi_job_is_array = False
        self.psi_job_id = None
        self.summary_job_id = None
        #self.psi_job_id = dict((sample_id, None) for sample_id in
        #                       self.sample_ids)
        #self.summary_job_id = dict((sample_id, None) for sample_id in
        #                           self.sample_ids)

        # self.psi_walltime = cl.args['psi_walltime']
        # self.summary_walltime = cl.args['summary_walltime']
        #
        # self.individual_jobs = cl.args['individual_jobs']

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

        #print 'self.read_type', self.read_type

        if self.read_type == 'paired_end':
            commands.append('\n# Calculate insert size')
            commands.append('python /projects/ps-yeolab/software/bin/pe_utils'
                            '.py '
                            '--compute-insert-len {0} '
                            '/projects/ps-yeolab/genomes/{1'
                            '}/miso/SE_constitutive/SE.{1}.min_20'
                            '.const_exons.gff '
                            '--no-bam-filter '
                            '--output-dir {2} '.format(bam, self.genome,
                                                       bam_dir))

            insert_len_stddev = 'INSERT_LEN_STDDEV'
            insert_len_mean = 'INSERT_LEN_MEAN'
            insert_len_file = bam + '.insert_len'
            commands.append(
                "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f1 | cut -d'=' -f2)"
                % (insert_len_mean, insert_len_file))

            # Assign {sample_id}_insert_len_STDDEV variable
            commands.append(
                "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f2 | cut -d'=' -f2)"
                % (insert_len_stddev, insert_len_file))

            insert_len_arguments = ' --paired-end $%s $%s ' % (insert_len_mean,
                                                               insert_len_stddev)

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