__author__ = 'olga'

from gscripts.qtools import Submitter
import os
from glob import glob
import sys
from gscripts import which
import pysam
from gscripts.general import read_sample_info_file


class MisoPipeline(object):
    def __init__(self, cl):
        """
        Given a CommandLine object 'cl', save the arguments as
        attributes of this class. And initialize the cluster job IDs as None,
         so we can check if they're there or not in the future.
        """
        self.read_type = cl.args['read_type']
        #self.event_type = cl.args['event_type'].upper()
        self.sample_info_file = cl.args['sample_info_file']

        self.debug = cl.args['debug']

        try:
            self.miso = (which('miso')[0])
            self.miso_scripts_dir = os.path.dirname(self.miso)
        except IndexError:
            # If there is an IndexError, that means that 'which' returned an
            # empty list, and thus there is no miso installed on the path.
            print >> sys.stderr, '"which miso" returned empty list, ' \
                                 'the program miso does not exist on your ' \
                                 'command line'
            sys.exit(1)

        self.paired_end_utils_py = '%s/pe_utils.py' % self.miso_scripts_dir

        # Remove the trailing slash. If it's not there, this won't do anything.
        self.annotation_index_strfmt = cl.args[
            'annotation_index_strfmt'].rstrip(
            '/')

        # Assuming we're using the annotation index structure described in
        # 'submit_miso_pipeline.py'
        # self.event_type_gff = glob(
        #     '%s/%s*.gff' % (self.base_annotation_dir, self.event_type))
        self.constitutive_exon_gff = cl.args['constitutive_exon_gff']
        #self.event_type_index = self.annotation_index_strfmt % self.event_type

        self.num_processes = cl.args['num_processes']
        self.num_cores = cl.args['num_cores']
        if self.sample_info_file:
            self.sample_ids, self.bams, self.notes = read_sample_info_file(
                self.sample_info_file)
        else:
            bam = cl.args['bam']
            sample_id = cl.args['sample_id']
            self.sample_ids = [sample_id]
            self.bams = [bam]

        self.output_sh = cl.args['output_sh']

        self.extra_miso_arguments = cl.args['extra_miso_arguments']

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

        self.psi_walltime = cl.args['psi_walltime']
        self.summary_walltime = cl.args['summary_walltime']

        self.individual_jobs = cl.args['individual_jobs']

        if cl.args['submit_sh_suffix'] != '':
            self.submit_sh_suffix = '_' + cl.args['submit_sh_suffix'].lstrip(
                '_')
        else:
            self.submit_sh_suffix = ''
        if cl.args['sample_id_suffix'] != '':
            self.sample_id_suffix = '_' + cl.args['sample_id_suffix'].lstrip(
                '_')
            self.sample_ids = [sample_id + self.sample_id_suffix
                               for sample_id in self.sample_ids]
        else:
            self.sample_id_suffix = ''

        self.sh_scripts_dir = cl.args['sh_scripts_dir'].rstrip('/')
        if self.sh_scripts_dir == '':
            self.sh_scripts_dir = os.curdir

        #self.job_name_prefix = 'miso%s_%s' % (self.submit_sh_suffix,
        #                                      self.event_type)

        self.queue = cl.args['queue']
        self.genome = cl.args['genome']


    def run_all_single_sample(self):

        bam = self.bams[0]
        sample_id = self.sample_ids[0]

        bam_dir = os.path.dirname(os.path.abspath(bam))

        commands = []
        commands.append('#!/bin/sh')
        commands.append('# Finding all MISO splicing scores for sample: {}. '
                        'Yay!\n'
                        .format(sample_id))

        insert_len_arguments = ''

        event_types = ['SE', 'MXE', 'AFE', 'ALE', 'A3SS', 'A5SS', 'TANDEMUTR']

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
            commands.append('python /home/yeo-lab/software/bin/pe_utils.py '
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
            commands.append('python /home/yeo-lab/software/bin/miso \
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


    def summary(self):
        summary_commands = []

        job_name_base = '%s_summary' % (self.job_name_prefix)
        submit_sh_base = '%s/%s' \
                         % (self.sh_scripts_dir, job_name_base)
        all_submit_sh = []

        for bam, sample_id, psi_output_dir, summary_output_dir in \
                zip(self.bams, self.sample_ids, self.psi_output_dirs,
                    self.summary_output_dirs):
            commands = []
            # Okay, now we are ready to write to the submitter script
            #commands.append('\n\n# --- %s --- #' % sample_id)

            # add a line of padding and the sample id to the output file
            #commands.append('\necho\necho "--- %s ----"' %
            #                sample_id)
            #commands.append('date')
            summary_command = 'python %s/run_miso.py --summarize-samples %s ' \
                              '%s >%s/summary.out 2>%s/summary.err' \
                              % (self.miso_scripts_dir, psi_output_dir,
                                 psi_output_dir, psi_output_dir,
                                 psi_output_dir)
            commands.append(summary_command)

            commands.append('# Copy over the summary files to prevent'
                            ' overloading the home directory')
            temp_summary_file = '%s/summary/%s.miso_summary' % (
                psi_output_dir, sample_id)
            final_summary_file = '%s/summary/%s.miso_summary' % (
                summary_output_dir, sample_id)
            commands.append('mkdir -p %s/summary' % (
                summary_output_dir.rstrip('/')))
            commands.append('cp %s %s' % (temp_summary_file,
                                          final_summary_file))

            # Join all these lines together into a single line, and add it to
            # the array job
            command = ' ; '.join(commands)
            summary_commands.append(command)
            if self.individual_jobs:
                job_name = '%s' % (job_name_base)
                submit_sh = '%s.sh' \
                            % (submit_sh_base)

                sh_out = submit_sh + '.out'
                sh_err = submit_sh + '.err'

                #if self.num_cores > 1:
                #    additional_resources = {'-t': '1-%d'
                #                                  % (self.num_processes *
                #                                     self.num_cores)}
                #else:
                #    additional_resources = None

                # if self.psi_job_id[sample_id] is not None:
                #     sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                #                     command_list=summary_commands,
                #                     job_name=job_name,
                #                     wait_for=[self.psi_job_id[sample_id]],
                #                     # Tell the queue to parallelize this job
                #                          # into a job array
                #                     additional_resources=additional_resources)
                # else:q
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                command_list=[command],
                                job_name=job_name,
                                # Tell the queue to parallelize this job
                                # into a job array. NOTE: changed the call of array
                                # job into where write_sh is called.
                                #additional_resources=additional_resources,
                                out=sh_out, err=sh_err)

                self.summary_job_id = sub.write_sh(
                    submit=True, nodes=self.num_cores, ppn=2, queue=self.queue,
                    walltime=self.summary_walltime)

                print self.summary_job_id

                # Put the submitter script wherever the command was run from
                #        if self.submit_sh_suffix:

                #        else:
                #            job_name = 'miso_%s_summary' % self.event_type
        job_name = '%s' % (job_name_base)
        submit_sh = '%s.sh' \
                    % (submit_sh_base)

        sh_out = submit_sh + '.out'
        sh_err = submit_sh + '.err'

        #if self.num_cores > 1:
        #    additional_resources = {'-t': '1-%d'
        #                                  % (self.num_processes *
        #                                     self.num_cores)}
        #else:
        #    additional_resources = None

        # if self.psi_job_id[sample_id] is not None:
        #     sub = Submitter(queue_type='PBS', sh_file=submit_sh,
        #                     command_list=summary_commands,
        #                     job_name=job_name,
        #                     wait_for=[self.psi_job_id[sample_id]],
        #                     # Tell the queue to parallelize this job
        #                          # into a job array
        #                     additional_resources=additional_resources)
        # else:
        sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                        command_list=summary_commands, job_name=job_name,
                        # Tell the queue to parallelize this job
                        # into a job array. NOTE: changed the call of array
                        # job into where write_sh is called.
                        #additional_resources=additional_resources,
                        out=sh_out, err=sh_err)

        self.summary_job_id = sub.write_sh(
            submit=True, nodes=self.num_cores, ppn=2, queue=self.queue,
            walltime=self.summary_walltime, array=True, max_running=20)

        print self.summary_job_id
        # Save all the qsub commands in one file
        #with open('%s.sh' % submit_sh_base, 'w') as f:
        #    # f.write('#!/bin/bash\n\n')
        #    f.writelines(all_submit_sh)